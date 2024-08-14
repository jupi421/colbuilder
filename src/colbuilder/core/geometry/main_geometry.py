from __future__ import annotations
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any, Union
import asyncio
import os
import subprocess
import shutil
from colorama import init, Fore, Style 

from colbuilder.core.geometry.crystal import Crystal
from colbuilder.core.geometry.crystalcontacts import CrystalContacts
from colbuilder.core.geometry.chimera import Chimera
from colbuilder.core.geometry.system import System
from colbuilder.core.geometry.connect import Connect
from colbuilder.core.geometry.caps import Caps
from colbuilder.core.geometry.optimize import Optimizer
from colbuilder.core.geometry.mix import Mix
from colbuilder.core.geometry.replace import Replace
from colbuilder.core.geometry.model import Model

from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import ColbuilderError
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

@timeit
async def build_geometry(config: ColbuilderConfig) -> Any:
    def ensure_pdb_extension(filename: str) -> str:
        if not filename.endswith('.pdb'):
            return filename + '.pdb'
        return filename

    try:
        path_wd = Path(config.working_directory)
        pdb_file = str(config.file)
        if pdb_file.endswith('.pdb'):
            pdb_file = pdb_file[:-4]
        contact_distance = config.contact_distance
        crystalcontacts_file = config.crystalcontacts_file
        connect_file = config.connect_file
        crystalcontacts_optimize = config.crystalcontacts_optimize
        solution_space = config.solution_space
        fibril_length = config.fibril_length
        geometry = config.geometry_generator
        pdb_out = Path(ensure_pdb_extension(str(config.output)))

        LOG.debug(f"Using PDB file: {pdb_file}")
        LOG.debug(f"Output PDB file: {pdb_out}")

        if pdb_file is None:
            raise ValueError('No PDB file given to build collagen microfibril')

        LOG.debug(f'Reading crystallographic symmetry from {pdb_file}')
        crystal = Crystal(pdb_file)
        crystal.translate_crystal(pdb=pdb_file, translate=[0, 0, 4000])

        path_pdb_file = path_wd / pdb_file
        chimera = Chimera(str(path_pdb_file))
        steps = 6

        if pdb_file and contact_distance and not crystalcontacts_file:
            LOG.info(f"Step 1/{steps} Building from contact distance")
            system, crystalcontacts, connect = build_from_contactdistance(
                path_wd, pdb_file, contact_distance, solution_space,
                crystalcontacts_file, chimera, crystal
            )
        elif pdb_file and not contact_distance and crystalcontacts_file:
            LOG.info(f"Step 1/{steps} Building from crystal contacts file")
            system, crystalcontacts, connect = build_from_crystalcontacts(
                crystalcontacts_file, solution_space, crystal, connect_file,
                crystalcontacts_optimize
            )
        elif not geometry:
            LOG.info('Set -geometry flag to generate microfibrillar structure PDB file')
            return system
        else:
            raise ValueError(
                'Please provide Contact Distance or CrystalContacts to generate the microfibril OR '
                'contact distance = 0 to generate the triple helix topology.'
            )

        if not geometry:
            LOG.info('Set -geometry flag to generate microfibrillar structure PDB file')
            return system

        LOG.debug(f'Writing {crystalcontacts.crystalcontacts_file}')
        crystalcontacts.write_crystalcontacts(system=system, crystalcontacts_file=crystalcontacts.crystalcontacts_file)

        LOG.info(f'Step 2/{steps} Generating system from {crystalcontacts.crystalcontacts_file}')
        LOG.info(f'{Fore.BLUE}Please wait, this may take some time ...{Style.RESET_ALL}')

        chimera.matrixset(pdb=pdb_file, crystalcontacts=crystalcontacts.crystalcontacts_file,
                          system_size=system.get_size(), fibril_length=fibril_length)

        LOG.info(f'Step 3/{steps} Cutting system to {fibril_length} nm')
        system = matrixset_system(system=system, crystalcontacts_file=crystalcontacts.crystalcontacts_file)

        LOG.info(f'Step 4/{steps} Writing {connect.connect_file}')
        connect.write_connect(system=system, connect_file=connect.connect_file)

        LOG.info(f'Step 5/{steps} Adding caps')
        rm_dir = os.path.join(path_wd, system.get_model(model_id=0.0).type)
        if os.path.exists(rm_dir):
            LOG.debug(f'Removing directory: {rm_dir}')
            try:
                shutil.rmtree(rm_dir)
            except Exception as e:
                LOG.warning(f'Failed to remove directory {rm_dir}: {str(e)}')
        else:
            LOG.debug(f'Directory does not exist, skipping removal: {rm_dir}')
        LOG.debug(f'Creating directory: {rm_dir}')
        try:
            os.makedirs(rm_dir, exist_ok=True)
        except Exception as e:
            LOG.error(f'Failed to create directory {rm_dir}: {str(e)}')
            raise ColbuilderError(f"Geometry generation failed: Unable to create directory. {rm_dir}")
        try:
            cap_system(system=system, crosslink_type=system.get_model(model_id=0.0).type)
        except Exception as e:
            LOG.error(f'Failed to add caps: {str(e)}')
            raise ColbuilderError(f"Geometry generation failed: Unable to add caps.")

        LOG.info(f'Step 6/{steps} Writing collagen fibril')
        system.write_pdb(pdb_out=Path(pdb_out), fibril_length=fibril_length)
        
        return system

    except FileNotFoundError as e:
        LOG.error(f"File not found error: {str(e)}")
        raise ColbuilderError(f"Geometry generation failed: {str(e)}")
    except Exception as e:
        LOG.error(f"An error occurred during geometry generation: {str(e)}")
        raise ColbuilderError(f"Geometry generation failed: {str(e)}")

@timeit
def build_system(crystal: Crystal, crystalcontacts: CrystalContacts) -> System:
    LOG.debug("Building system of models")
    system = System(crystal=crystal, crystalcontacts=crystalcontacts)

    transformation = system.crystalcontacts.read_t_matrix()
    unit_cell: Dict[float, Any] = {
        k: system.crystal.get_s_matrix(t_matrix=transformation[k]) 
        for k in transformation
    }

    for key_m in transformation:
        model = Model(
            id=key_m,
            transformation=transformation[key_m],
            unit_cell=unit_cell[key_m],
            pdb_file=crystal.pdb_file
        )
        system.add_model(model=model)

    LOG.debug(f"Built system with {len(system.get_models())} models")
    return system

@timeit
def build_from_contactdistance(
    path_wd: Path, 
    pdb_file: str, 
    contact_distance: float,
    solution_space: List[float], 
    crystalcontacts_file: Optional[str],
    chimera: Chimera, 
    crystal: Crystal
) -> Tuple[System, CrystalContacts, Connect]:
    path_pdb_file = path_wd / pdb_file
    crystalcontacts_file = 'crystalcontacts_from_colbuilder'  # default name
    connect_file = 'connect_from_colbuilder'  # default name
    
    LOG.info(f'     Getting CrystalContacts for contact distance {contact_distance} Ang')
    chimera.matrixget(pdb=str(path_pdb_file),
                      contact_distance=contact_distance,
                      crystalcontacts=crystalcontacts_file)
    
    LOG.info(f'     Writing {crystalcontacts_file}')
    crystalcontacts = CrystalContacts(crystalcontacts_file)
    
    LOG.info(f'     Building system')
    system = build_system(crystal=crystal, crystalcontacts=crystalcontacts)
    
    LOG.info(f'     Connecting system')
    system, connect = connect_system(system=system, connect_file=connect_file)
    
    LOG.info(f'     Optimizing system')
    optimizer = Optimizer(system=system, solution_space=solution_space)
    system = optimizer.run_optimize(system=system, connect=connect)
    system, connect = connect_system(system=system, connect_file=connect_file)
    
    crystalcontacts.crystalcontacts_file = crystalcontacts_file + '_opt'
    
    return system, crystalcontacts, connect

@timeit
def build_from_crystalcontacts(
    crystalcontacts_file: Union[str, Path],
    solution_space: List[float],
    crystal: Crystal, 
    connect_file: Optional[Union[str, Path]],
    crystalcontacts_optimize: bool
) -> Tuple[System, CrystalContacts, Connect]:
    crystalcontacts_file = Path(crystalcontacts_file).with_suffix('')
    if not connect_file:
        bool_external_connect = False
        connect_file = 'connect_from_colbuilder'  # default name
    else:
        bool_external_connect = True
        connect_file = Path(connect_file).with_suffix('')
        
    crystalcontacts = CrystalContacts(str(crystalcontacts_file))
    
    LOG.info(f'     Building system')
    system = build_system(crystal=crystal, crystalcontacts=crystalcontacts)
    
    LOG.info(f'     Connecting system')
    system, connect = connect_system(system=system, connect_file=connect_file)
    
    if crystalcontacts_optimize:
        LOG.info(f'     Optimizing system')
        optimizer = Optimizer(system=system, solution_space=solution_space)
        system = optimizer.run_optimize(system=system, connect=connect)
        system, connect = connect_system(system=system)
        
        crystalcontacts.crystalcontacts_file = crystalcontacts.crystalcontacts_file + '_opt'
    
    if bool_external_connect:
        system = connect.get_external_connect_file(system=system, connect_file=connect_file)
        
    return system, crystalcontacts, connect

@timeit
def connect_system(system: System, connect_file: Optional[str] = None) -> Tuple[System, Connect]:
    LOG.debug("Identifying crosslink connections")
    connect = Connect(system=system, connect_file=connect_file)
    system_connect = connect.run_connect(system=system)
    for key_m in system_connect:
        system.get_model(model_id=key_m).add_connect(
                connect_id=key_m,
                connect=system_connect[key_m]
            )
    LOG.debug("Crosslink connections identified")
    return system, connect  

@timeit
def cap_system(system: System, crosslink_type: str) -> Caps:
    LOG.debug("Capping models in the system with terminal groups")
    caps = Caps(system=system)
    for idx in system.get_models():
        pdb_id = int(idx)
        pdb_file = f"{pdb_id}.pdb"
        if not os.path.exists(pdb_file):
            LOG.warning(f"PDB file {pdb_file} not found. Model IDs: {list(system.get_models())}")
            continue
        caps.read_residues(pdb_id=pdb_id)
        caps.add_caps(pdb_id=pdb_id, crosslink_type=crosslink_type)
    LOG.debug("Caps added to models")
    return caps

@timeit
def matrixset_system(system: System, crystalcontacts_file: Union[str, Path]) -> System:
    LOG.debug(f"Setting system after cutting fibril. crystalcontacts_file type: {type(crystalcontacts_file)}")
    crystalcontacts_file = Path(crystalcontacts_file)
    LOG.debug(f"crystalcontacts_file after conversion: {crystalcontacts_file}")
    if crystalcontacts_file.suffix == '.txt':
        crystalcontacts_file = crystalcontacts_file.with_suffix('') 
    id_file = crystalcontacts_file.with_name(f"{crystalcontacts_file.name}_id.txt")
    LOG.debug(f"id_file: {id_file}")
     
    if not id_file.exists():
        raise FileNotFoundError(f"Crystal contacts ID file not found: {id_file}")
    
    try:
        with open(id_file, 'r') as f:
            contacts = [float(i.split(' ')[1]) for i in f.readlines()]
    except Exception as e:
        LOG.error(f"Error reading crystal contacts ID file: {str(e)}")
        raise

    for model in list(system.get_models()):  # Create a list to avoid modifying during iteration
        if model not in contacts:
            system.delete_model(model_id=model)
        elif system.get_model(model_id=model).connect is not None:
            for connect in list(system.get_model(model_id=model).connect):  # Create a list to avoid modifying during iteration
                if connect not in contacts: 
                    system.get_model(model_id=model).delete_connect(connect_id=connect)
    
    LOG.debug("System set after cutting")    
    return system

@timeit
async def mix_geometry(system: System, config: ColbuilderConfig) -> System:
    path_wd = Path(config.working_directory)
    fibril_length = config.fibril_length
    connect_file = config.connect_file
    pdb_files = [Path(pdb) for pdb in config.files_mix]
    ratio_mix = config.ratio_mix
    pdb_out = Path(config.output)
    steps = 2

    if not pdb_files:
        raise ValueError("No PDB files provided for mixing")

    try:
        if ratio_mix and not connect_file:
            mix_setup = ratio_mix
            mix_pdb = dict(zip(mix_setup.keys(), pdb_files))
            
            system_size = system.get_size()
            connect_file = Path('connect_from_colbuilder')
            
            LOG.info(f'Step 1/{steps} Generating mix setup')
            for key in list(mix_pdb.keys()):
                Crystal(pdb=str(mix_pdb[key])).translate_crystal(pdb=str(mix_pdb[key]), translate=[0, 0, 4000])
                chimera_ = Chimera(str(path_wd / mix_pdb[key]))
                
                LOG.info(f' System {key}:')
                LOG.info(f'     Generating system from {mix_pdb[key]} {Fore.BLUE}...{Style.RESET_ALL}')
                chimera_.matrixset(pdb=str(mix_pdb[key]), crystalcontacts=str(system.crystalcontacts.crystalcontacts_file),
                                system_size=system_size, fibril_length=fibril_length)

                LOG.info(f'     Cutting system to {fibril_length} nm')
                system = matrixset_system(system=system, crystalcontacts_file=str(system.crystalcontacts.crystalcontacts_file))
                
                LOG.info(f'     Adding caps')
                dir_path = path_wd / key
                if dir_path.exists():
                    shutil.rmtree(dir_path)
                dir_path.mkdir(exist_ok=True)
                cap_system(system=system, crosslink_type=key)
            
            LOG.info(f'Step 2/{steps} Mixing systems')
            mix_ = Mix(ratio_mix=mix_setup, system=system)
            system = mix_.add_mix(system=system)
            
        elif not ratio_mix and connect_file:
            for pdb in pdb_files[1:]:
                Crystal(pdb=str(pdb)).translate_crystal(pdb=str(pdb), translate=[0, 0, 4000])
            
            connect_file = Path(connect_file) if connect_file else None
            
            LOG.info(f'{Fore.YELLOW}NOTE: Make sure each crosslink-type in {connect_file} (D,T,DT,TD) is provided{Style.RESET_ALL}')
            mix_ = Mix(system=system, connect_mix=str(connect_file) if connect_file else None)
            system = mix_.get_mix_from_connect_file(connect_file=str(connect_file) if connect_file else None)
        else:
            raise ValueError('Either provide a connect-file with the crosslink-type OR '
                             'provide the ratio_mix and files_mix to generate a differently crosslinked microfibril')

        connect_file_path = Path(connect_file) if connect_file else Path('connect_from_colbuilder.txt')
        LOG.debug(f'Writing connect file {connect_file_path}')
        Connect(system=system).write_connect(system=system, connect_file=connect_file_path)
        
        pdb_out_path = Path(pdb_out) if pdb_out else Path(f'{pdb_out}.pdb')
        system.write_pdb(pdb_out=pdb_out_path, fibril_length=fibril_length)
        
    except Exception as e:
        LOG.error(f"Error in mixing process: {str(e)}")
        LOG.error(f"Error type: {type(e)}")
        import traceback
        LOG.error(f"Traceback: {traceback.format_exc()}")
        raise ValueError(f"Error in mixing process: {str(e)}")

    LOG.info(f'{Fore.BLUE}Mixing geometry process completed{Style.RESET_ALL}')
    return system

@timeit
async def replace_geometry(system: System, config: ColbuilderConfig) -> System:
    path_wd = Path(config.working_directory)
    ratio_replace = config.ratio_replace
    fibril_length = config.fibril_length
    pdb_out = Path(config.output)
    replace_file = config.replace_file
    steps = 2

    try:
        if not replace_file:
            replace = Replace(ratio_replace=ratio_replace, system=system, fibril_length=fibril_length)
            system, current_ratio = replace.run_replace(system=system, ratio_replace=ratio_replace)
            LOG.info(f"Step 1/{steps} Writing replace file")
            replace.write_replace(system=system, file='replace')
            replace_file = 'replace'
            LOG.info(f'Step 2/{steps} Replacing {current_ratio:.4f}% of crosslinks from collagen fibril by LYS')
        else:
            LOG.info(f'Step 1/{steps} Using existing replace file: {replace_file}')
            LOG.info(f'Step 2/{steps} Replacing {ratio_replace}% of crosslinks from collagen fibril by LYS')

        LOG.info(f'{Fore.BLUE}Please wait, this may take some time ...{Style.RESET_ALL}')
        
        chimera = Chimera(str(path_wd / system.crystal.pdb_file))
        
        result = chimera.swapaa(replace=replace_file, system_type=system.get_model(model_id=0.0).type)
        
        if result.returncode != 0:
            LOG.error(f"Chimera swapaa command failed with return code {result.returncode}")
            raise RuntimeError(f"Chimera swapaa command failed. Check logs for details.")

        pdb_out_path = Path(pdb_out) if pdb_out else Path(f'{pdb_out}.pdb')
        system.write_pdb(pdb_out=pdb_out_path, fibril_length=fibril_length)

    except Exception as e:
        LOG.error(f"Error in replace_geometry: {str(e)}")
        import traceback
        LOG.error(f"Traceback: {traceback.format_exc()}")
        raise ValueError(f"Error in replace_geometry: {str(e)}")

    LOG.info(f'{Fore.BLUE}Replace geometry process completed.{Style.RESET_ALL}')
    return system

def main():
    pass

if __name__ == "__main__":
    main()
