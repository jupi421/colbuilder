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
        
    async def replace(self, config: ColbuilderConfig) -> System:
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
            path_wd = Path(config.working_directory)
            ratio_replace = config.ratio_replace
            fibril_length = config.fibril_length
            pdb_out = Path(config.output)
            replace_file = config.replace_file
            steps = 2
            LOG.info(f'Step 1/{steps} Using existing replace file: {replace_file}')
            LOG.info(f'Step 2/{steps} Replacing {ratio_replace}% of crosslinks from collagen fibril by LYS')
            
            LOG.info(f'{Fore.BLUE}Please wait, this may take some time ...{Style.RESET_ALL}')
            
            chimera = Chimera(config, str(path_wd / system.crystal.pdb_file))
        
            result = chimera.swapaa(replace=replace_file, system_type=system.get_model(model_id=0.0).type)
        
            if result.returncode != 0:
                LOG.error(f"Chimera swapaa command failed with return code {result.returncode}")
                raise RuntimeError(f"Chimera swapaa command failed. Check logs for details.")

            pdb_out_path = Path(pdb_out) if pdb_out else Path(f'{pdb_out}.pdb')
            system.write_pdb(pdb_out=pdb_out_path, fibril_length=fibril_length)
                    
        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Failed to replace geometry elements",
                original_error=e,
                error_code="GEO_ERR_004",
                context={"replace_file": str(config.replace_file)}
            )
            
