from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Optional, Union, Any, Tuple
import shutil
from colorama import Fore, Style
import os

from colbuilder.core.geometry.crystal import Crystal
from colbuilder.core.geometry.system import System
from colbuilder.core.geometry.chimera import Chimera
from colbuilder.core.geometry.connect import Connect
from colbuilder.core.geometry.caps import Caps
from colbuilder.core.geometry.crystalcontacts import CrystalContacts
from colbuilder.core.geometry.mix import Mix
from colbuilder.core.geometry.optimize import Optimizer
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.exceptions import ColbuilderError

LOG = setup_logger(__name__)

class CrosslinkMixer:
    """Class for mixing different crosslink types in collagen systems."""
    
    def __init__(self):
        LOG.info("info: Initializing CrosslinkMixer")
        self.path_wd: Optional[Path] = None
        self.fibril_length: Optional[float] = None
        self.contact_distance: Optional[float] = None
        self.crystalcontacts_file: str = 'crystalcontacts_from_colbuilder'
        LOG.info("info: CrosslinkMixer initialized successfully")
    
    @staticmethod
    def _ensure_pdb_extension(filename: str) -> str:
        """Ensure filename has .pdb extension."""
        if not filename.endswith('.pdb'):
            return filename + '.pdb'
        return filename
    
    @staticmethod
    def _build_system(crystal: Crystal, crystalcontacts: CrystalContacts) -> System:
        """Build a system from crystal and crystal contacts."""
        LOG.info("Building system of models")
        system = System(crystal=crystal, crystalcontacts=crystalcontacts)
        
        transformation = system.crystalcontacts.read_t_matrix()
        unit_cell: Dict[float, Any] = {
            k: system.crystal.get_s_matrix(t_matrix=transformation[k])
            for k in transformation
        }
        
        from colbuilder.core.geometry.model import Model
        for key_m in transformation:
            model = Model(
                id=key_m,
                transformation=transformation[key_m],
                unit_cell=unit_cell[key_m],
                pdb_file=crystal.pdb_file
            )
            system.add_model(model=model)
        
        LOG.info(f"Built system with {len(system.get_models())} models")
        return system
    
    def _build_from_contactdistance(
        self,
        path_wd: Path,
        pdb_file: str,
        contact_distance: float,
        solution_space: List[float],
        crystalcontacts_file: Optional[str],
        chimera: Chimera,
        crystal: Crystal
    ) -> Tuple[System, CrystalContacts, Connect]:
        """Build system from contact distance."""
        path_pdb_file = path_wd / pdb_file
        connect_file = 'connect_from_colbuilder'
        
        LOG.info(f'     Getting CrystalContacts for contact distance {contact_distance} Ang')
        chimera.matrixget(
            pdb=str(path_pdb_file),
            contact_distance=contact_distance,
            crystalcontacts=self.crystalcontacts_file
        )
        
        LOG.info(f'     Writing {self.crystalcontacts_file}')
        crystalcontacts = CrystalContacts(self.crystalcontacts_file)
        
        LOG.info(f'     Building system')
        system = self._build_system(crystal=crystal, crystalcontacts=crystalcontacts)
        
        LOG.info(f'     Connecting system')
        system, connect = self._connect_system(system=system, connect_file=connect_file)
        
        has_crosslinks = any(
            hasattr(system.get_model(model_id=model_id), 'crosslink') and
            system.get_model(model_id=model_id).crosslink
            for model_id in system.get_models()
        )
        
        if has_crosslinks:
            LOG.info(f'     Optimizing system')
            optimizer = Optimizer(system=system, solution_space=solution_space)
            system = optimizer.run_optimize(system=system, connect=connect)
            system, connect = self._connect_system(system=system, connect_file=connect_file)
            crystalcontacts.crystalcontacts_file = self.crystalcontacts_file + '_opt'
        else:
            LOG.info(f'     Skipping optimization for non-crosslinked system')
        
        return system, crystalcontacts, connect
    
    @staticmethod
    def _connect_system(system: System, connect_file: Optional[str] = None) -> Tuple[System, Connect]:
        """Connect models in the system."""
        LOG.info("Identifying crosslink connections")
        connect = Connect(system=system, connect_file=connect_file)
        system_connect = connect.run_connect(system=system)
        
        for key_m in system_connect:
            system.get_model(model_id=key_m).add_connect(
                connect_id=key_m,
                connect=system_connect[key_m]
            )
        
        LOG.info("Crosslink connections identified")
        return system, connect
    
    @staticmethod
    def _cap_system(system: System, crosslink_type: str) -> Caps:
        """Add caps to the system."""
        LOG.info("Capping models in the system with terminal groups")
        caps = Caps(system=system)
        
        for idx in system.get_models():
            pdb_id = int(idx)
            pdb_file = f"{pdb_id}.pdb"
            if not os.path.exists(pdb_file):
                LOG.warning(f"PDB file {pdb_file} not found. Model IDs: {list(system.get_models())}")
                continue
            caps.read_residues(pdb_id=pdb_id)
            caps.add_caps(pdb_id=pdb_id, crosslink_type=crosslink_type)
            
        LOG.info("Caps added to models")
        return caps
    
    def _matrixset_system(self, system: System) -> System:
        """Update system after cutting to specified length."""
        LOG.info(f"Setting system after cutting fibril")
        
        id_file = Path(f"{self.crystalcontacts_file}_opt_id.txt")
        LOG.info(f"Looking for ID file: {id_file}")
        
        if not id_file.exists():
            raise FileNotFoundError(f"Crystal contacts ID file not found: {id_file}")
        
        try:
            with open(id_file, 'r') as f:
                contacts = [float(i.split(' ')[1]) for i in f.readlines()]
        except Exception as e:
            LOG.error(f"Error reading crystal contacts ID file: {str(e)}")
            raise
        
        initial_models = len(system.get_models())
        for model in list(system.get_models()):
            if model not in contacts:
                system.delete_model(model_id=model)
            elif system.get_model(model_id=model).connect is not None:
                for connect in list(system.get_model(model_id=model).connect):
                    if connect not in contacts:
                        system.get_model(model_id=model).delete_connect(connect_id=connect)
        
        LOG.info(f"System cut from {initial_models} to {len(system.get_models())} models")
        return system

    async def mix(self, system: Optional[System], config: ColbuilderConfig) -> System:
        """Main method to mix different crosslink types."""
        LOG.info("info: Starting mix method")
        try:
            self.path_wd = Path(config.working_directory)
            self.fibril_length = config.fibril_length
            self.contact_distance = config.contact_distance
            
            # Set up mixing parameters first
            mix_setup = config.ratio_mix
            mix_pdb = dict(zip(mix_setup.keys(), config.files_mix))
            
            LOG.info(f'Step 1/2 Generating mix setup')
            
            # Check if we need to initialize system
            if system.get_size() == 0:
                first_pdb = str(config.files_mix[0])
                if first_pdb.endswith('.pdb'):
                    first_pdb = first_pdb[:-4]
                
                LOG.info(f"No system provided, building initial system from {first_pdb}")
                crystal = Crystal(first_pdb)
                crystal.translate_crystal(pdb=first_pdb, translate=[0, 0, 4000])
                
                chimera = Chimera(config, str(self.path_wd / first_pdb))
                
                # Build initial system
                system, crystalcontacts, connect = self._build_from_contactdistance(
                    self.path_wd,
                    first_pdb,
                    self.contact_distance,
                    config.solution_space,
                    None,
                    chimera,
                    crystal
                )
                
                # Write crystalcontacts file
                LOG.info(f'Writing {crystalcontacts.crystalcontacts_file}')
                crystalcontacts.write_crystalcontacts(
                    system=system, 
                    crystalcontacts_file=crystalcontacts.crystalcontacts_file
                )
                
                LOG.info(f"Initial system built, proceeding with mixing")
            else:
                LOG.info("Using provided system")
            
            system_size = system.get_size()
            connect_file = Path('connect_from_colbuilder')
            
            LOG.info(f"System size: {system_size}")
            LOG.info(f"Crystalcontacts file: {system.crystalcontacts.crystalcontacts_file}")
            
            for key in list(mix_setup.keys()):
                Crystal(pdb=str(mix_pdb[key])).translate_crystal(
                    pdb=str(mix_pdb[key]),
                    translate=[0, 0, 4000]
                )
                chimera = Chimera(config, str(self.path_wd / mix_pdb[key]))
                
                LOG.info(f' System {key}:')
                LOG.info(f'     Generating system from {mix_pdb[key]} {Fore.BLUE}...{Style.RESET_ALL}')
                
                chimera.matrixset(
                    pdb=str(mix_pdb[key]),
                    crystalcontacts=str(system.crystalcontacts.crystalcontacts_file),
                    system_size=system_size,
                    fibril_length=self.fibril_length
                )
                
                LOG.info(f'     Cutting system to {self.fibril_length} nm')
                system = self._matrixset_system(system)
                
                LOG.info(f'     Adding caps')
                dir_path = self.path_wd / key
                if dir_path.exists():
                    shutil.rmtree(dir_path)
                dir_path.mkdir(exist_ok=True)
                self._cap_system(system, key)
            
            LOG.info(f'Step 2/2 Mixing systems')
            mix_ = Mix(ratio_mix=mix_setup, system=system)
            system = mix_.add_mix(system=system)
            
            connect_file_path = Path('connect_from_colbuilder.txt')
            LOG.info(f'Writing connect file {connect_file_path}')
            Connect(system=system).write_connect(system=system, connect_file=connect_file_path)
        
            pdb_out_path = Path(config.output) if config.output else Path(f'{config.output}.pdb')
            system.write_pdb(pdb_out=pdb_out_path, fibril_length=self.fibril_length)
        
            return system
            
        except Exception as e:
            LOG.info(f"info: Exception in mix: {str(e)}")
            LOG.info(f"info: Exception type: {type(e)}")
            try:
                os.chdir(self.path_wd)
            except:
                pass
            LOG.info("Full state at error:")
            LOG.info(f"Has system: {system is not None}")
            if system:
                LOG.info(f"System size: {system.get_size()}")
                LOG.info(f"Has crystalcontacts: {hasattr(system, 'crystalcontacts')}")
            raise