# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any
import asyncio
import os
import subprocess
import shutil
from colorama import init, Fore, Style

from colbuilder.geometry.crystal import Crystal
from colbuilder.geometry.crystalcontacts import CrystalContacts
from colbuilder.geometry.chimera import Chimera
from colbuilder.geometry.system import System
from colbuilder.geometry.connect import Connect
from colbuilder.geometry.caps import Caps
from colbuilder.geometry.optimize import Optimizer
from colbuilder.geometry.mix import Mix
from colbuilder.geometry.fibril import Fibril
from colbuilder.geometry.replace import Replace
from colbuilder.geometry.model import Model

from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import ColbuilderError
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

@timeit
async def build_geometry(config: ColbuilderConfig) -> Any:
    """
    Build system of models from input configuration.

    Parameters
    ----------
    config : ColbuilderConfig
        Configuration object containing all necessary parameters.

    Returns
    -------
    System
        Built system of models.

    Raises
    ------
    ValueError
        If no PDB file is provided or if invalid configuration is given.
    """
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
        elif pdb_file and contact_distance == 0 and not crystalcontacts_file:
            LOG.info(f"Step 1/{steps} Building from pdb file")
            system, crystalcontacts, connect = build_from_pdb(
                path_wd, pdb_file, contact_distance,
                crystalcontacts_file, chimera, crystal
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
            raise ColbuilderError(f"Geometry generation failed: Unable to create directory {rm_dir}")
        try:
            cap_system(system=system, crosslink_type=system.get_model(model_id=0.0).type)
        except Exception as e:
            LOG.error(f'Failed to add caps: {str(e)}')
            raise ColbuilderError(f"Geometry generation failed: Unable to add caps")

        LOG.info(f'Step 6/{steps} Writing collagen fibril.')
        system.write_pdb(pdb_out=Path(pdb_out), fibril_length=fibril_length)

        system.write_pdb(pdb_out=Path(pdb_out), fibril_length=fibril_length)
        
        return system

    except FileNotFoundError as e:
        LOG.error(f"File not found error: {str(e)}")
        raise ColbuilderError(f"Geometry generation failed: {str(e)}")
    except Exception as e:
        LOG.error(f"An error occurred during geometry generation: {str(e)}")
        raise ColbuilderError(f"Geometry generation failed: {str(e)}")

@timeit
def replace_geometry(config: Dict[str, Any], system: System) -> System:
    """
    Replace crosslinks within microfibril to reduce the overall number of crosslinks.

    Parameters
    ----------
    config : Dict[str, Any]
        Configuration dictionary.
    system : System
        The system to modify.

    Returns
    -------
    System
        Modified system with replaced crosslinks.
    """
    path_wd = Path(config.working_directory)
    ratio_replace = config.ratio_replace
    fibril_length = config.fibril_length
    pdb_out = config.output
    replace_file = config.replace_file

    if not replace_file:
        LOG.info(f'Replacing {ratio_replace}% of crosslinks from microfibril')
        replace = Replace(ratio_replace=ratio_replace, system=system, fibril_length=fibril_length)
        system = replace.run_replace(system=system, ratio_replace=ratio_replace)

        replace.write_replace(system=system, file='replace')
        replace_file = 'replace'
    else:
        system = system

    LOG.info('Please wait, this may take some time ...')
    chimera = Chimera(str(path_wd / system.crystal.pdb_file))
    chimera.swapaa(replace=replace_file, system_type=system.get_model(model_id=0.0).type)

    system.write_pdb(pdb_out=pdb_out, fibril_length=fibril_length)

    return system

@timeit
def mix_geometry(config: Dict[str, Any], system: System) -> System:
    """
    Mix differently crosslinked models/triple helices within the system.

    Parameters
    ----------
    config : Dict[str, Any]
        Configuration dictionary.
    system : System
        The system to modify.

    Returns
    -------
    System
        Modified system with mixed crosslinks.

    Raises
    ------
    ValueError
        If invalid configuration is provided.
    """
    path_wd = Path(config.working_directory)
    fibril_length = config.fibril_length
    connect_file = config.connect_file
    pdb_files = config.files_mix
    ratio_mix = config.ratio_mix
    pdb_out = config.output

    LOG.info('Preparing mix setup')
    if ratio_mix and not connect_file:
        mix_setup = {idx.split(':')[0]: idx.split(':')[1] for idx in ratio_mix}
        mix_pdb = dict(zip(mix_setup.keys(), pdb_files))

        system_size = system.get_size()
        connect_file = 'connect_from_colbuilder'

        for key in list(mix_pdb.keys())[1:]:
            Crystal(pdb=mix_pdb[key]).translate_crystal(pdb=mix_pdb[key], translate=[0, 0, 4000])
            chimera = Chimera(str(path_wd / mix_pdb[key]))

            LOG.info(f'Generating {key} system from {mix_pdb[key]}')
            LOG.info('Please wait, this may take some time ...')
            chimera.matrixset(pdb=mix_pdb[key], crystalcontacts=system.crystalcontacts.crystalcontacts_file,
                              system_size=system_size, fibril_length=fibril_length)

            LOG.info(f'Cutting {key} system to {fibril_length} nm')
            system = matrixset_system(system=system, crystalcontacts_file=system.crystalcontacts.crystalcontacts_file)

            LOG.info('Adding caps')
            subprocess.run(f'rm -r {path_wd}/{key}', shell=True)
            subprocess.run(f'mkdir {path_wd}/{key}', shell=True)
            cap_system(system=system, crosslink_type=key)

        LOG.info('Mixing system')
        mix = Mix(ratio_mix=mix_setup, system=system)
        system = mix.add_mix(system=system)

    elif not ratio_mix and connect_file:
        for pdb in pdb_files[1:]:
            Crystal(pdb=pdb).translate_crystal(pdb=pdb, translate=[0, 0, 4000])
        connect_file = connect_file.replace('.txt', '')

        LOG.info('Mixing system')
        LOG.info(f'NOTE: Make sure each crosslink-type in {connect_file} (D,T,DT,TD) is provided')
        mix = Mix(system=system, connect_mix=connect_file)
        system = mix.get_mix_from_connect_file(connect_file=connect_file)
    else:
        raise ValueError('Either provide a connect-file with the crosslink-type OR '
                         'provide the -ratio_mix and -files_mix to generate a differently crosslinked microfibril')

    Connect(system=system).write_connect(system=system, connect_file=connect_file)

    system.write_pdb(pdb_out=pdb_out, fibril_length=fibril_length)

    return system

@timeit
def build_fibril(path_wd: str, pdb_file: str, connect_file: Optional[str]) -> System:
    """
    Build a system from colbuilder 1.0 fibril.

    Parameters
    ----------
    path_wd : str
        Working directory path.
    pdb_file : str
        PDB file name.
    connect_file : Optional[str]
        Connect file name.

    Returns
    -------
    System
        Built system from the fibril.
    """
    system = System(pdb_fibril=pdb_file)

    LOG.info(f'Reading fibril {pdb_file} from colbuilder 1.0')
    fibril = Fibril(system=system, pdb_file=pdb_file)

    LOG.info('Separating system')
    fibril.separate_system(pdb_file=pdb_file)

    LOG.info('Building system')
    system = fibril.build_system(system=system)

    LOG.info('Connecting system')
    system, _ = connect_system(system=system)

    LOG.info('Writing connect_from_colbuilder.txt')
    fibril.write_connect(system=system, connect_file='connect_from_colbuilder')

    return system

@timeit
def build_system(crystal: Crystal, crystalcontacts: CrystalContacts) -> System:
    """
    Build a system of models, i.e., a microfibril of triple helices.

    Parameters
    ----------
    crystal : Crystal
        The crystal object containing crystal information.
    crystalcontacts : CrystalContacts
        The crystal contacts object.

    Returns
    -------
    System
        A system of models representing the microfibril.
    """
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
def connect_system(system: System, connect_file: Optional[str] = None) -> Tuple[System, Connect]:
    """
    Identify crosslink connections between models within the system.

    Parameters
    ----------
    system : System
        The system to connect.
    connect_file : Optional[str], default=None
        Path to the connect file, if any.

    Returns
    -------
    Tuple[System, Connect]
        Connected system and the Connect object.
    """
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
    """
    Cap each model of system.
    
    Parameters
    ----------
    system : System
        The system to cap.
    crosslink_type : str
        Type of crosslink to use for capping.
        
    Returns
    -------
    Caps
        The Caps object used for capping.
    """
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
def matrixset_system(system: System, crystalcontacts_file: str) -> System:
    """
    Set system after cutting the fibril to its desired length.
    
    Parameters
    ----------
    system : System
        The system to modify.
    crystalcontacts_file : str
        Path to the crystal contacts file.
        
    Returns
    -------
    System
        The modified system.
    """
    LOG.debug("Setting system after cutting fibril")
    id_file = crystalcontacts_file + '_id.txt'
    if not os.path.exists(id_file):
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
def build_from_contactdistance(
    path_wd: Path, 
    pdb_file: str, 
    contact_distance: float,
    solution_space: List[float], 
    crystalcontacts_file: Optional[str],
    chimera: Chimera, 
    crystal: Crystal
) -> Tuple[System, CrystalContacts, Connect]:
    """
    Generate system of models based on contact distance and PDB-file.
    
    Parameters
    ----------
    path_wd : Path
        Working directory path.
    pdb_file : str
        PDB file name.
    contact_distance : float
        Contact distance for crystal contacts.
    solution_space : List[float]
        Solution space parameters.
    crystalcontacts_file : Optional[str]
        Crystal contacts file name, if any.
    chimera : Chimera
        Chimera object for operations.
    crystal : Crystal
        Crystal object.
        
    Returns
    -------
    Tuple[System, CrystalContacts, Connect]
        Built system, crystal contacts, and connections.
    """
    path_pdb_file = path_wd / pdb_file
    crystalcontacts_file = 'crystalcontacts_from_colbuilder' # default name
    connect_file = 'connect_from_colbuilder' # default name
    
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
    
    crystalcontacts.crystalcontacts_file = crystalcontacts_file+'_opt'
    
    return system, crystalcontacts, connect

@timeit
def build_from_crystalcontacts(
    crystalcontacts_file: str, 
    solution_space: List[float],
    crystal: Crystal, 
    connect_file: Optional[str],
    crystalcontacts_optimize: bool
) -> Tuple[System, CrystalContacts, Connect]:
    """
    Generate system of models based on CrystalContacts and PDB-file.
    
    Parameters
    ----------
    crystalcontacts_file : str
        Crystal contacts file name.
    solution_space : List[float]
        Solution space parameters.
    crystal : Crystal
        Crystal object.
    connect_file : Optional[str]
        Connect file name, if any.
    crystalcontacts_optimize : bool
        Whether to optimize crystal contacts.
        
    Returns
    -------
    Tuple[System, CrystalContacts, Connect]
        Built system, crystal contacts, and connections.
    """
    LOG.info(f"Building from crystal contacts")
    crystalcontacts_file = crystalcontacts_file.replace('.txt','')
    if not connect_file:
        bool_external_connect = False
        connect_file = 'connect_from_colbuilder' # default name
    else:
        bool_external_connect = True
        connect_file = connect_file.replace('.txt','')
        
    crystalcontacts = CrystalContacts(crystalcontacts_file)
    steps = 3
    
    LOG.info(f'Step 1/{steps} Building system')
    system = build_system(crystal=crystal, crystalcontacts=crystalcontacts)
    
    LOG.info(f'Step 2/{steps}Connecting system')
    system, connect = connect_system(system=system, connect_file=connect_file)
    
    if crystalcontacts_optimize:
        LOG.info(f'Step 3/{steps}Optimizing system')
        optimizer = Optimizer(system=system, solution_space=solution_space)
        system = optimizer.run_optimize(system=system, connect=connect)
        system, connect = connect_system(system=system)
        
        crystalcontacts.crystalcontacts_file = crystalcontacts.crystalcontacts_file+'_opt'
    
    if bool_external_connect:
        system = connect.get_external_connect_file(system=system, connect_file=connect_file)
        
    LOG.debug("Build from crystal contacts completed")
    return system, crystalcontacts, connect

@timeit
def build_from_pdb(
    path_wd: Path, 
    pdb_file: str, 
    contact_distance: float,
    crystalcontacts_file: Optional[str], 
    chimera: Chimera,
    crystal: Crystal
) -> Tuple[System, CrystalContacts, Connect]:
    """
    Generate system of models based on contact distance and PDB-file.
    
    Parameters
    ----------
    path_wd : Path
        Working directory path.
    pdb_file : str
        PDB file name.
    contact_distance : float
        Contact distance for crystal contacts.
    crystalcontacts_file : Optional[str]
        Crystal contacts file name, if any.
    chimera : Chimera
        Chimera object for operations.
    crystal : Crystal
        Crystal object.
        
    Returns
    -------
    Tuple[System, CrystalContacts, Connect]
        Built system, crystal contacts, and connections.
    """
    LOG.info(f"{Fore.BLUE}Building from PDB{Style.RESET_ALL}")
    
    path_pdb_file = path_wd / pdb_file
    
    crystalcontacts_file = 'crystalcontacts_from_colbuilder' # default name
    steps = 3
    
    LOG.info(f'Step 1/{steps} Preparing topology for single PDB-file')
    chimera.matrixget(pdb=str(path_pdb_file), contact_distance=contact_distance,
                      crystalcontacts=crystalcontacts_file)
    
    crystalcontacts = CrystalContacts(crystalcontacts_file)
    
    LOG.info(f'Step 2/{steps} Building system')
    system = build_system(crystal=crystal, crystalcontacts=crystalcontacts)
    system, connect = connect_system(system=system)
    
    LOG.debug(f"Step 3/{steps} Build from PDB completed")
    return system, crystalcontacts, connect

def main():
    pass

if __name__ == "__main__":
    main()