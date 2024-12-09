"""
Colbuilder Mixer Module

This module provides services for mixing different geometry types in the system.
It supports mixing based on ratios or connect files.

Key Features:
    - Ratio-based mixing
    - Connect file-based mixing
    - System validation
    - PDB file processing
    - Resource management
"""

from pathlib import Path
from typing import Dict, List, Optional, Any
import shutil
import os
from colorama import Fore, Style

from ..utils.exceptions import GeometryGenerationError
from ..utils.logger import setup_logger
from ..utils.config import ColbuilderConfig

from .crystal import Crystal
from .system import System
from .mix import Mix
from .chimera import Chimera
from .caps import Caps

LOG = setup_logger(__name__)

class CrosslinkMixer:
    """
    Service for mixing geometry types in the system.
    
    This class manages the process of combining different geometries
    based on mixing ratios or connection files.
    """
    
    def __init__(self):
        """Initialize the mixer service."""
        self.steps = 2

    async def mix(self, system: System, config: ColbuilderConfig) -> System:
        """
        Mix different geometry types in the system.
        
        Args:
            system: System to mix
            config: Configuration settings
            
        Returns:
            System: Mixed system
            
        Raises:
            GeometryGenerationError: If mixing fails
        """
        try:
            if config.ratio_mix and not config.connect_file:
                return await self._mix_with_ratio(system, config)
            elif not config.ratio_mix and config.connect_file:
                return await self._mix_with_connect_file(system, config)
            else:
                raise GeometryGenerationError(
                    message="Invalid mixing configuration",
                    error_code="GEO_ERR_003",
                    context={
                        "ratio_mix": str(config.ratio_mix),
                        "connect_file": str(config.connect_file)
                        if config.connect_file else None
                    }
                )
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to mix system",
                original_error=e,
                error_code="GEO_ERR_003",
                context={"config": config.model_dump()}
            )
            
    async def _mix_with_ratio(self, system: System, config: ColbuilderConfig) -> System:
        """
        Mix system using ratio configuration.
        
        Args:
            system: System to mix
            config: Configuration settings
            
        Returns:
            System: Mixed system
            
        Raises:
            GeometryGenerationError: If ratio mixing fails
        """
        try:
            mix_setup = config.ratio_mix
            pdb_files = [Path(pdb) for pdb in config.files_mix]
            
            if not pdb_files:
                raise GeometryGenerationError(
                    message="No PDB files provided for mixing",
                    error_code="GEO_ERR_003",
                    context={"files_mix": str(config.files_mix)}
                )
            
            mix_pdb = dict(zip(mix_setup.keys(), pdb_files))
            system_size = system.get_size()
            connect_file = Path('connect_from_colbuilder')
            
            LOG.info(f'Step 1/{self.steps} Generating mix setup')
            await self._process_mix_components(mix_pdb, system, system_size, config)
            
            LOG.info(f'Step 2/{self.steps} Mixing systems')
            try:
                mixer = Mix(ratio_mix=mix_setup, system=system)
                system = mixer.add_mix(system=system)
            except Exception as e:
                raise GeometryGenerationError(
                    message="Failed to mix system components",
                    original_error=e,
                    error_code="GEO_ERR_003",
                    context={"ratio_mix": str(mix_setup)}
                )
            
            await self._finalize_mixing(system, connect_file, config)
            return system
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to mix system with ratios",
                original_error=e,
                error_code="GEO_ERR_003",
                context={
                    "ratio_mix": str(config.ratio_mix),
                    "pdb_files": [str(p) for p in pdb_files]
                }
            )

    async def _process_mix_components(
        self,
        mix_pdb: Dict[str, Path],
        system: System,
        system_size: Any,
        config: ColbuilderConfig
    ) -> None:
        """
        Process individual mixing components.
        
        Args:
            mix_pdb: Dictionary mapping mix types to PDB files
            system: Current system
            system_size: Size of the system
            config: Configuration settings
            
        Raises:
            GeometryGenerationError: If component processing fails
        """
        for mix_type, pdb_path in mix_pdb.items():
            try:
                await self._process_single_component(
                    mix_type,
                    pdb_path,
                    system,
                    system_size,
                    config
                )
            except Exception as e:
                raise GeometryGenerationError(
                    message=f"Failed to process mix component {mix_type}",
                    original_error=e,
                    error_code="GEO_ERR_003",
                    context={
                        "mix_type": mix_type,
                        "pdb_path": str(pdb_path)
                    }
                )
                
    async def _process_single_component(
        self,
        mix_type: str,
        pdb_path: Path,
        system: System,
        system_size: Any,
        config: ColbuilderConfig
    ) -> None:
        """
        Process a single mixing component.
        
        Args:
            mix_type: Type of mix component
            pdb_path: Path to PDB file
            system: Current system
            system_size: Size of the system
            config: Configuration settings
            
        Raises:
            GeometryGenerationError: If processing fails
        """
        try:
            crystal = Crystal(pdb=str(pdb_path))
            crystal.translate_crystal(
                pdb=str(pdb_path),
                translate=[0, 0, 4000]
            )
            
            chimera = Chimera(
                config,
                str(Path(config.working_directory) / pdb_path)
            )
            
            LOG.info(f' System {mix_type}:')
            LOG.info(
                f'     Generating system from {pdb_path} '
                f'{Fore.BLUE}...{Style.RESET_ALL}'
            )
            
            await self._generate_component_matrix(
                chimera,
                pdb_path,
                system,
                system_size,
                config
            )
            
            await self._add_component_caps(mix_type, system, config)
            
        except Exception as e:
            raise GeometryGenerationError(
                message=f"Failed to process component {mix_type}",
                original_error=e,
                error_code="GEO_ERR_003",
                context={
                    "mix_type": mix_type,
                    "pdb_path": str(pdb_path),
                    "system_size": system_size
                }
            )

    async def _generate_component_matrix(
        self,
        chimera: Chimera,
        pdb_path: Path,
        system: System,
        system_size: Any,
        config: ColbuilderConfig
    ) -> None:
        """
        Generate matrix for a mixing component.
        
        Args:
            chimera: Chimera instance
            pdb_path: Path to PDB file
            system: Current system
            system_size: Size of the system
            config: Configuration settings
            
        Raises:
            GeometryGenerationError: If matrix generation fails
        """
        try:
            chimera.matrixset(
                pdb=str(pdb_path),
                crystalcontacts=str(system.crystalcontacts.crystalcontacts_file),
                system_size=system_size,
                fibril_length=config.fibril_length
            )
            
            LOG.info(f"     Cutting system to {config.fibril_length} nm")
            await self._cut_component_system(
                system,
                system.crystalcontacts.crystalcontacts_file
            )
            
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to generate component matrix",
                original_error=e,
                error_code="GEO_ERR_003",
                context={
                    "pdb_path": str(pdb_path),
                    "system_size": system_size,
                    "fibril_length": config.fibril_length
                }
            )

    async def _cut_component_system(
        self,
        system: System,
        crystalcontacts_file: str
    ) -> System:
        """
        Cut component system to specified length.
        
        Args:
            system: System to cut
            crystalcontacts_file: Crystal contacts file path
            
        Returns:
            System: Cut system
            
        Raises:
            GeometryGenerationError: If cutting fails
        """
        try:
            crystalcontacts_path = Path(crystalcontacts_file)
            
            if crystalcontacts_path.suffix == '.txt':
                crystalcontacts_path = crystalcontacts_path.with_suffix('')
                
            id_file = crystalcontacts_path.with_name(
                f"{crystalcontacts_path.name}_id.txt"
            )
            
            if not id_file.exists():
                raise GeometryGenerationError(
                    message="Crystal contacts ID file not found",
                    error_code="GEO_ERR_003",
                    context={"id_file": str(id_file)}
                )
            
            try:
                with open(id_file, 'r') as f:
                    contacts = [float(i.split(' ')[1]) for i in f.readlines()]
            except Exception as e:
                raise GeometryGenerationError(
                    message="Failed to read crystal contacts ID file",
                    original_error=e,
                    error_code="GEO_ERR_003",
                    context={"id_file": str(id_file)}
                )

            for model_id in list(system.get_models()):
                if model_id not in contacts:
                    system.delete_model(model_id=model_id)
                elif system.get_model(model_id=model_id).connect is not None:
                    for connect_id in list(system.get_model(model_id=model_id).connect):
                        if connect_id not in contacts:
                            system.get_model(model_id=model_id).delete_connect(
                                connect_id=connect_id
                            )
            
            return system
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to cut component system",
                original_error=e,
                error_code="GEO_ERR_003",
                context={"crystalcontacts_file": str(crystalcontacts_file)}
            )

    async def _add_component_caps(
        self,
        mix_type: str,
        system: System,
        config: ColbuilderConfig
    ) -> None:
        """
        Add caps to component system.
        
        Args:
            mix_type: Type of mix component
            system: System to cap
            config: Configuration settings
            
        Raises:
            GeometryGenerationError: If capping fails
        """
        try:
            LOG.info("     Adding caps")
            dir_path = Path(config.working_directory) / mix_type
            
            if dir_path.exists():
                shutil.rmtree(dir_path)
            dir_path.mkdir(exist_ok=True)
            
            caps = Caps(system=system)
            caps.add_caps(system=system, crosslink_type=mix_type)
            
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to add component caps",
                original_error=e,
                error_code="GEO_ERR_003",
                context={
                    "mix_type": mix_type,
                    "dir_path": str(dir_path)
                }
            )
            
    async def _mix_with_connect_file(
        self,
        system: System,
        config: ColbuilderConfig
    ) -> System:
        """
        Mix system using connect file.
        
        Args:
            system: System to mix
            config: Configuration settings
            
        Returns:
            System: Mixed system
            
        Raises:
            GeometryGenerationError: If connect file mixing fails
        """
        try:
            pdb_files = [Path(pdb) for pdb in config.files_mix]
            
            for pdb_path in pdb_files[1:]:
                try:
                    Crystal(pdb=str(pdb_path)).translate_crystal(
                        pdb=str(pdb_path),
                        translate=[0, 0, 4000]
                    )
                except Exception as e:
                    raise GeometryGenerationError(
                        message=f"Failed to translate crystal for {pdb_path}",
                        original_error=e,
                        error_code="GEO_ERR_003",
                        context={"pdb_path": str(pdb_path)}
                    )
                    
            connect_file = (
                Path(config.connect_file) if config.connect_file else None
            )
            
            LOG.info(
                f'{Fore.YELLOW}NOTE: Make sure each crosslink-type in '
                f'{connect_file} (D,T,DT,TD) is provided{Style.RESET_ALL}'
            )
            
            try:
                mixer = Mix(
                    system=system,
                    connect_mix=str(connect_file) if connect_file else None
                )
                system = mixer.get_mix_from_connect_file(
                    connect_file=str(connect_file) if connect_file else None
                )
            except Exception as e:
                raise GeometryGenerationError(
                    message="Failed to mix using connect file",
                    original_error=e,
                    error_code="GEO_ERR_003",
                    context={"connect_file": str(connect_file)}
                )
            
            await self._finalize_mixing(system, connect_file, config)
            return system
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to mix system with connect file",
                original_error=e,
                error_code="GEO_ERR_003",
                context={
                    "connect_file": str(config.connect_file),
                    "pdb_files": [str(p) for p in pdb_files]
                }
            )
            
    async def _finalize_mixing(
        self,
        system: System,
        connect_file: Optional[Path],
        config: ColbuilderConfig
    ) -> None:
        """
        Finalize the mixing process.
        
        Args:
            system: Mixed system
            connect_file: Optional connect file path
            config: Configuration settings
            
        Raises:
            GeometryGenerationError: If finalization fails
        """
        try:
            connect_file_path = (
                Path(connect_file) if connect_file
                else Path('connect_from_colbuilder.txt')
            )
            
            try:
                LOG.debug(f'Writing connect file {connect_file_path}')
                from .connect import Connect
                Connect(system=system).write_connect(
                    system=system,
                    connect_file=connect_file_path
                )
            except Exception as e:
                raise GeometryGenerationError(
                    message="Failed to write connect file",
                    original_error=e,
                    error_code="GEO_ERR_003",
                    context={"connect_file": str(connect_file_path)}
                )
            
            try:
                pdb_out_path = (
                    Path(config.output) if config.output
                    else Path(f'{config.output}.pdb')
                )
                system.write_pdb(
                    pdb_out=pdb_out_path,
                    fibril_length=config.fibril_length
                )
            except Exception as e:
                raise GeometryGenerationError(
                    message="Failed to write output PDB file",
                    original_error=e,
                    error_code="GEO_ERR_003",
                    context={
                        "output_file": str(pdb_out_path),
                        "fibril_length": config.fibril_length
                    }
                )
            
            LOG.info(f'{Fore.BLUE}Mixing geometry process completed{Style.RESET_ALL}')
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to finalize mixing",
                original_error=e,
                error_code="GEO_ERR_003",
                context={"connect_file": str(connect_file) if connect_file else None}
            )