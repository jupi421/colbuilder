"""
Colbuilder Geometry Replacer Module

This module implements the GeometryReplacer operation, which handles the
replacement of specific elements in the geometry system, such as replacing
crosslinks with lysines.

The replacer supports two modes:
1. Replacement using a specified ratio
2. Replacement using an existing replace file

Key Features:
    - Ratio-based replacement
    - File-based replacement
    - Amino acid swapping
    - System validation
    - Progress tracking
    - Resource management

Example Usage:
    replacer = GeometryReplacer(context)
    modified_system = await replacer.execute(system)
"""

from pathlib import Path
from typing import Optional, Tuple
from colorama import Fore, Style

from ..utils.exceptions import GeometryGenerationError
from ..utils.logger import setup_logger
from ..utils.config import ColbuilderConfig

from .system import System
from .replace import Replace
from .chimera import Chimera

LOG = setup_logger(__name__)

class GeometryReplacer:
    """
    Handles replacement of geometry elements in the system.
    
    This class manages the process of replacing specific elements
    in the geometry, such as replacing crosslinks with lysines.
    It supports both ratio-based and file-based replacement approaches.
    
    Attributes:
        context: Shared geometry context
        steps: Number of steps in the replacement process
    """
    
    def __init__(self):
        """Initialize the geometry replacer."""
        self.steps = 2
        
    async def replace(self, system: System, config: ColbuilderConfig) -> System:
        """
        Replace crosslink residues by their standard amino acid versions.
        
        Args:
            system: Current system state
            
        Returns:
            System: System with replacements
            
        Raises:
            GeometryGenerationError: If replacement fails or no system provided
        """
        try:
            replace_file = await self._prepare_replacement(system, config)
            return await self._perform_replacement(system, replace_file, config)
                    
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to replace geometry elements",
                original_error=e,
                error_code="GEO_ERR_004",
                context={"replace_file": str(config.replace_file)}
            )
            
    async def _prepare_replacement(self, system: System, config: ColbuilderConfig) -> Path:
        """
        Prepare system for replacement operation.
        
        Args:
            system: System to prepare
            
        Returns:
            Path: Path to replacement file
            
        Raises:
            GeometryGenerationError: If preparation fails
        """
        try:
            if not config.replace_file:
                return await self._prepare_ratio_replacement(system, config)
            else:
                LOG.info(
                    f'Step 1/{self.steps} Using existing replace file: '
                    f'{config.replace_file}'
                )
                LOG.info(
                    f'Step 2/{self.steps} Replacing {config.ratio_replace}% '
                    'of crosslinks from collagen fibril by LYS'
                )
                return Path(config.replace_file)
                
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to prepare replacement",
                original_error=e,
                error_code="GEO_ERR_004",
                context={
                    "ratio_replace": config.ratio_replace,
                    "replace_file": str(config.replace_file)
                    if config.replace_file else None
                }
            )
            
    async def _prepare_ratio_replacement(
        self,
        system: System,
        config: ColbuilderConfig
    ) -> Tuple[Path, float]:
        """
        Prepare ratio-based replacement.
        
        Args:
            system: System to prepare
            
        Returns:
            Path: Path to generated replace file
            
        Raises:
            GeometryGenerationError: If preparation fails
        """
        try:
            replacer = Replace(
                ratio_replace=config.ratio_replace,
                system=system,
                fibril_length=config.fibril_length
            )
            
            try:
                system, current_ratio = replacer.run_replace(
                    system=system,
                    ratio_replace=config.ratio_replace
                )
            except Exception as e:
                raise GeometryGenerationError(
                    message="Failed to calculate replacement ratio",
                    original_error=e,
                    error_code="GEO_ERR_004",
                    context={
                        "ratio_replace": config.ratio_replace,
                        "fibril_length": config.fibril_length
                    }
                )
            
            LOG.info(f"Step 1/{self.steps} Writing replace file")
            try:
                replacer.write_replace(system=system, file='replace')
            except Exception as e:
                raise GeometryGenerationError(
                    message="Failed to write replace file",
                    original_error=e,
                    error_code="GEO_ERR_004",
                    context={"file": 'replace'}
                )
            
            LOG.info(
                f'Step 2/{self.steps} Replacing {current_ratio:.4f}% of '
                'crosslinks from collagen fibril by LYS'
            )
            
            return Path('replace')
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to prepare ratio replacement",
                original_error=e,
                error_code="GEO_ERR_004",
                context={"ratio_replace": config.ratio_replace}
            )
            
    async def _perform_replacement(
        self,
        system: System,
        replace_file: Path, 
        config: ColbuilderConfig
    ) -> System:
        """
        Perform the actual replacement operation.
        
        Args:
            system: System to modify
            replace_file: Path to replacement file
            
        Returns:
            System: Modified system
            
        Raises:
            GeometryGenerationError: If replacement fails
        """
        try:
            LOG.info(f'{Fore.BLUE}Please wait, this may take some time ...{Style.RESET_ALL}')
            
            chimera = await self._initialize_chimera(system, config)
            
            await self._swap_amino_acids(chimera, replace_file, system)
            
            await self._write_final_structure(system, config)
            
            LOG.info(f'{Fore.BLUE}Replace geometry process completed.{Style.RESET_ALL}')
            return system
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to perform replacement",
                original_error=e,
                error_code="GEO_ERR_004",
                context={
                    "pdb_file": str(system.crystal.pdb_file),
                    "replace_file": str(replace_file)
                }
            )
            
    async def _initialize_chimera(self, system: System, config: ColbuilderConfig) -> Chimera:
        """
        Initialize Chimera for replacement.
        
        Args:
            system: Current system
            
        Returns:
            Chimera: Initialized Chimera instance
            
        Raises:
            GeometryGenerationError: If initialization fails
        """
        try:
            return Chimera(
                config,
                str(Path(config.working_directory) / system.crystal.pdb_file)
            )
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to initialize Chimera",
                original_error=e,
                error_code="GEO_ERR_004",
                context={"pdb_file": str(system.crystal.pdb_file)}
            )
            
    async def _swap_amino_acids(
        self,
        chimera: Chimera,
        replace_file: Path,
        system: System
    ) -> None:
        """
        Perform amino acid swap operation.
        
        Args:
            chimera: Initialized Chimera instance
            replace_file: Path to replacement file
            system: Current system
            
        Raises:
            GeometryGenerationError: If swap operation fails
        """
        try:
            result = chimera.swapaa(
                replace=replace_file,
                system_type=system.get_model(model_id=0.0).type
            )
            
            if result.returncode != 0:
                raise GeometryGenerationError(
                    message="Chimera swap amino acids command failed",
                    error_code="GEO_ERR_004",
                    context={
                        "return_code": result.returncode,
                        "replace_file": str(replace_file),
                        "system_type": system.get_model(model_id=0.0).type
                    }
                )
                
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to swap amino acids",
                original_error=e,
                error_code="GEO_ERR_004",
                context={
                    "replace_file": str(replace_file),
                    "system_type": system.get_model(model_id=0.0).type
                }
            )
            
    async def _write_final_structure(self, system: System, config: ColbuilderConfig) -> None:
        """
        Write the final system structure.
        
        Args:
            system: Modified system
            
        Raises:
            GeometryGenerationError: If writing fails
        """
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
                message="Failed to write final structure",
                original_error=e,
                error_code="GEO_ERR_004",
                context={
                    "output_file": str(pdb_out_path),
                    "fibril_length": config.fibril_length
                }
            )