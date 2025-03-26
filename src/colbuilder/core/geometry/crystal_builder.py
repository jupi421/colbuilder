"""
Colbuilder Crystal Builder Module

This module provides services for building and managing crystal systems,
including initialization, contact distance calculations, and system optimization.

Key Features:
    - Crystal structure initialization
    - Contact distance calculations
    - System optimization
    - Cap addition
    - PDB file generation
"""

from pathlib import Path
from typing import Tuple, Optional, List, Dict, Any, Union
import shutil
import os
from colorama import Fore, Style

from ..utils.exceptions import GeometryGenerationError
from ..utils.logger import setup_logger
from ..utils.config import ColbuilderConfig
from ..utils.dec import timeit
from ..utils.files import FileManager

from .crystal import Crystal
from .system import System
from .model import Model
from .crystalcontacts import CrystalContacts
from .connect import Connect
from .caps import Caps
from .optimize import Optimizer
from .chimera import Chimera

LOG = setup_logger(__name__)

class CrystalBuilder:
    """
    Service for crystal system operations.
    
    This class encapsulates all operations related to building and managing
    crystal systems, including initialization, optimization, and file handling.
    """
    
    def __init__(self, file_manager: Optional[FileManager] = None):
        """Initialize the crystal builder."""
        self.steps = 7
        self.file_manager = file_manager
    
    def set_file_manager(self, file_manager: FileManager) -> None:
        """Set the file manager for consistent file handling."""
        self.file_manager = file_manager

    @timeit
    def matrixset_system(self, system: System, crystalcontacts_file: Union[str, Path]) -> System:
        """
        Process and update system after matrix generation.
        
        Args:
            system: System to process
            crystalcontacts_file: Path to crystal contacts file
            
        Returns:
            System: Processed system
            
        Raises:
            FileNotFoundError: If required files are missing
            GeometryGenerationError: If processing fails
        """
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

        for model in list(system.get_models()):  
            if model not in contacts:
                system.delete_model(model_id=model)
            elif system.get_model(model_id=model).connect is not None:
                for connect in list(system.get_model(model_id=model).connect):  
                    if connect not in contacts: 
                        system.get_model(model_id=model).delete_connect(connect_id=connect)
        
        LOG.debug("System set after cutting")    
        return system
        
    async def build(self, config: ColbuilderConfig) -> System:
        """Build a crystal system from configuration."""
        if not self.file_manager:
            self.file_manager = FileManager(config)
            
        try:
            LOG.info(f"Step 1/{self.steps} Initializing crystal")
            crystal = await self._initialize_crystal(config)
            
            LOG.info(f"Step 2/{self.steps} Building initial system")
            system, crystalcontacts, connect = await self._build_initial_system(crystal, config)
            
            LOG.info(f"Step 3/{self.steps} Writing crystal contacts")
            crystalcontacts.write_crystalcontacts(
                system=system,
                crystalcontacts_file=crystalcontacts.crystalcontacts_file
            )
            
            LOG.info(f"Step 4/{self.steps} Generating system matrix")
            LOG.info(f'{Fore.BLUE}Please wait, this may take some time ...{Style.RESET_ALL}')
            
            pdb_path = Path(config.pdb_file).resolve() if config.pdb_file else None
            chimera_scripts_dir = self.file_manager.find_file('chimera_scripts')
            
            chimera = Chimera(config, str(pdb_path))
            chimera.matrixset(
                pdb=str(config.pdb_file),
                crystalcontacts=crystalcontacts.crystalcontacts_file,
                system_size=system.get_size(),
                fibril_length=config.fibril_length
            )
            
            system = self.matrixset_system(system=system, crystalcontacts_file=crystalcontacts.crystalcontacts_file)
           
            LOG.info(f"Step 5/{self.steps} Writing {connect.connect_file}")
            connect.write_connect(system=system, connect_file=connect.connect_file)
            
            LOG.info(f"Step 6/{self.steps} Adding caps")
            model_type = system.get_model(model_id=0.0).type
            
            type_dir = self.file_manager.get_type_dir(model_type)
            
            has_crosslinks = self._needs_optimization(system)
            if has_crosslinks:
                caps = Caps(system=system)
                for idx in system.get_models():
                    pdb_id = int(idx)
                    pdb_file = f"{pdb_id}.pdb"
                    if not os.path.exists(pdb_file):
                        LOG.warning(f"PDB file {pdb_file} not found. Model IDs: {list(system.get_models())}")
                        continue
                    caps.read_residues(pdb_id=pdb_id)
                    caps.add_caps(pdb_id=pdb_id, crosslink_type=model_type)
            else:
                LOG.debug("System has no crosslinks, using standard capping")
                caps = Caps(system=system)
                for idx in system.get_models():
                    pdb_id = int(idx)
                    pdb_file = f"{pdb_id}.pdb"
                    if not os.path.exists(pdb_file):
                        LOG.warning(f"PDB file {pdb_file} not found. Model IDs: {list(system.get_models())}")
                        continue
                    caps.read_residues(pdb_id=pdb_id)
                    caps.add_caps(pdb_id=pdb_id, crosslink_type="NC")
            
            LOG.info(f"Step 7/{self.steps} Writing final structure")
            output_path = self.file_manager.get_output_path(config.output, ".pdb")
            system.write_pdb(
                pdb_out=output_path,
                fibril_length=config.fibril_length,
                cleanup=False
            )
            
            return system
            
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to build crystal system",
                original_error=e,
                error_code="GEO_ERR_001",
                context={"config": config.model_dump()}
            )
            
    async def _initialize_crystal(self, config: ColbuilderConfig) -> Crystal:
        """
        Initialize crystal structure from PDB file.
        
        Args:
            config: Configuration settings
            
        Returns:
            Crystal: Initialized crystal
            
        Raises:
            GeometryGenerationError: If initialization fails
        """
        try:
            if not config.pdb_file:
                raise GeometryGenerationError(
                    message="PDB file not specified",
                    error_code="GEO_ERR_005",
                    context={"config": config.model_dump()}
                )
                
            pdb_path = Path(config.pdb_file).resolve()
            if not pdb_path.exists():
                raise GeometryGenerationError(
                    message=f"PDB file not found: {pdb_path}",
                    error_code="GEO_ERR_005"
                )
                
            pdb_str = str(pdb_path)
            if pdb_str.endswith('.pdb'):
                pdb_str = pdb_str[:-4]
                
            original_pdb = self.file_manager.get_output_path(f"{os.path.basename(pdb_str)}_original", ".pdb")
            try:
                shutil.copy(f"{pdb_str}.pdb", original_pdb)
            except (shutil.Error, IOError) as e:
                raise GeometryGenerationError(
                    message="Failed to create backup of original PDB file",
                    original_error=e,
                    error_code="GEO_ERR_005",
                    context={
                        "source": f"{pdb_str}.pdb",
                        "destination": str(original_pdb)
                    }
                )
                
            crystal = Crystal(pdb_str)
            crystal.translate_crystal(
                pdb=pdb_str,
                translate=[0, 0, 4000]
            )
            return crystal
            
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to initialize crystal structure",
                original_error=e,
                error_code="GEO_ERR_002",
                context={"pdb_file": str(pdb_path)}
            )
    
    async def _build_initial_system(
        self,
        crystal: Crystal,
        config: ColbuilderConfig
    ) -> Tuple[System, CrystalContacts, Connect]:
        """
        Build initial system from crystal structure.
        
        Args:
            crystal: Crystal structure
            config: Configuration settings
            
        Returns:
            Tuple containing:
                - System: Built system
                - CrystalContacts: Crystal contacts information
                - Connect: System connections
                
        Raises:
            GeometryGenerationError: If system building fails
        """
        if config.contact_distance and not config.crystalcontacts_file:
            return await self._build_from_contact_distance(crystal, config)
        elif not config.contact_distance and config.crystalcontacts_file:
            return await self._build_from_crystal_contacts(crystal, config)
        else:
            raise GeometryGenerationError(
                message="Invalid configuration for system building",
                error_code="GEO_ERR_001",
                context={
                    "contact_distance": config.contact_distance,
                    "crystalcontacts_file": str(config.crystalcontacts_file)
                    if config.crystalcontacts_file else None
                }
            )

    async def _build_from_contact_distance(
        self,
        crystal: Crystal,
        config: ColbuilderConfig
    ) -> Tuple[System, CrystalContacts, Connect]:
        """
        Build system using contact distance method.
        """
        try:
            pdb_path = Path(config.working_directory) / config.pdb_file
            
            chimera = Chimera(config, str(pdb_path))
            crystalcontacts_file = 'crystalcontacts_from_colbuilder'
            connect_file = 'connect_from_colbuilder'
            
            LOG.debug(f'Getting CrystalContacts for contact distance {config.contact_distance} Ang')
            
            chimera.matrixget(
                pdb=str(pdb_path),
                contact_distance=config.contact_distance,
                crystalcontacts=crystalcontacts_file
            )
            
            crystalcontacts = CrystalContacts(crystalcontacts_file)
            system = await self._build_system_structure(crystal, crystalcontacts)
            system, connect = await self._connect_system(system, connect_file)
            
            if self._needs_optimization(system):
                LOG.info(f'{Fore.BLUE}Optimizing system{Style.RESET_ALL}')
                system = await self._optimize_system(
                    system,
                    connect,
                    config.solution_space
                )
                crystalcontacts.crystalcontacts_file = crystalcontacts_file + '_opt'
            
            return system, crystalcontacts, connect
            
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to build system from contact distance",
                original_error=e,
                error_code="GEO_ERR_002",
                context={
                    "contact_distance": config.contact_distance,
                    "pdb_file": str(pdb_path)
                }
            )

    async def _build_from_crystal_contacts(
        self,
        crystal: Crystal,
        config: ColbuilderConfig
    ) -> Tuple[System, CrystalContacts, Connect]:
        """
        Build system from existing crystal contacts file.
        
        Args:
            crystal: Crystal structure
            config: Configuration settings
            
        Returns:
            Tuple containing system components
            
        Raises:
            GeometryGenerationError: If building fails
        """
        try:
            crystalcontacts_file = Path(config.crystalcontacts_file).with_suffix('')
            connect_file = (Path(config.connect_file) if config.connect_file 
                          else Path('connect_from_colbuilder'))
            
            try:
                crystalcontacts = CrystalContacts(str(crystalcontacts_file))
            except Exception as e:
                raise GeometryGenerationError(
                    message="Failed to load crystal contacts file",
                    original_error=e,
                    error_code="GEO_ERR_002",
                    context={"crystalcontacts_file": str(crystalcontacts_file)}
                )
            
            system = await self._build_system_structure(crystal, crystalcontacts)
            system, connect = await self._connect_system(system, str(connect_file))
            
            if config.crystalcontacts_optimize:
                LOG.info(f'{Fore.BLUE}Optimizing system{Style.RESET_ALL}')
                system = await self._optimize_system(
                    system,
                    connect,
                    config.solution_space
                )
                crystalcontacts.crystalcontacts_file = str(crystalcontacts_file) + '_opt'
            
            if config.connect_file:
                try:
                    system = connect.get_external_connect_file(
                        system=system,
                        connect_file=str(connect_file)
                    )
                except Exception as e:
                    raise GeometryGenerationError(
                        message="Failed to process external connect file",
                        original_error=e,
                        error_code="GEO_ERR_002",
                        context={"connect_file": str(connect_file)}
                    )
            
            return system, crystalcontacts, connect
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to build system from crystal contacts",
                original_error=e,
                error_code="GEO_ERR_002",
                context={
                    "crystalcontacts_file": str(config.crystalcontacts_file),
                    "connect_file": str(config.connect_file) if config.connect_file else None
                }
            )

    async def _build_system_structure(
        self,
        crystal: Crystal,
        crystalcontacts: CrystalContacts
    ) -> System:
        """
        Build core system structure from crystal and contacts.
        
        Args:
            crystal: Crystal structure
            crystalcontacts: Crystal contacts information
            
        Returns:
            System: Built system
            
        Raises:
            GeometryGenerationError: If structure building fails
        """
        try:
            system = System(crystal=crystal, crystalcontacts=crystalcontacts)
            transformation = system.crystalcontacts.read_t_matrix()
            
            unit_cell = {
                model_id: system.crystal.get_s_matrix(t_matrix=transformation[model_id])
                for model_id in transformation
            }
            
            for model_id in transformation:
                try:
                    model = Model(
                        id=model_id,
                        transformation=transformation[model_id],
                        unit_cell=unit_cell[model_id],
                        pdb_file=crystal.pdb_file
                    )
                    system.add_model(model=model)
                except Exception as e:
                    raise GeometryGenerationError(
                        message="Failed to create and add model",
                        original_error=e,
                        error_code="GEO_ERR_002",
                        context={"model_id": model_id}
                    )
            
            return system
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to build system structure",
                original_error=e,
                error_code="GEO_ERR_002",
                context={
                    "crystal_pdb": crystal.pdb_file,
                    "crystalcontacts_file": crystalcontacts.crystalcontacts_file
                }
            )

    async def _connect_system(
        self,
        system: System,
        connect_file: str
    ) -> Tuple[System, Connect]:
        """
        Establish connections in the system.
        
        Args:
            system: System to connect
            connect_file: Connection file path
            
        Returns:
            Tuple containing connected system and connection information
            
        Raises:
            GeometryGenerationError: If connection fails
        """
        try:
            LOG.debug("Identifying crosslink connections")
            connect = Connect(system=system, connect_file=connect_file)
            system_connect = connect.run_connect(system=system)
            
            for model_id in system_connect:
                try:
                    system.get_model(model_id=model_id).add_connect(
                        connect_id=model_id,
                        connect=system_connect[model_id]
                    )
                except Exception as e:
                    raise GeometryGenerationError(
                        message="Failed to add connection to model",
                        original_error=e,
                        error_code="GEO_ERR_002",
                        context={"model_id": model_id}
                    )
            
            return system, connect
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to connect system",
                original_error=e,
                error_code="GEO_ERR_002",
                context={"connect_file": connect_file}
            )

    def _needs_optimization(self, system: System) -> bool:
        """
        Check if system needs optimization.
        
        Args:
            system: System to check
            
        Returns:
            bool: True if system has crosslinks
        """
        return any(
            hasattr(system.get_model(model_id), 'crosslink') and 
            system.get_model(model_id).crosslink 
            for model_id in system.get_models()
        )

    async def _optimize_system(
        self,
        system: System,
        connect: Connect,
        solution_space: List[float]
    ) -> System:
        """
        Optimize system structure.
        
        Args:
            system: System to optimize
            connect: System connections
            solution_space: Optimization parameters
            
        Returns:
            System: Optimized system
            
        Raises:
            GeometryGenerationError: If optimization fails
        """
        try:
            optimizer = Optimizer(system=system, solution_space=solution_space)
            optimized_system = optimizer.run_optimize(system=system, connect=connect)
            optimized_system, _ = await self._connect_system(
                optimized_system,
                connect.connect_file
            )
            return optimized_system
            
        except Exception as e:
            raise GeometryGenerationError(
                message="System optimization failed",
                original_error=e,
                error_code="GEO_ERR_002",
                context={"solution_space": solution_space}
            )