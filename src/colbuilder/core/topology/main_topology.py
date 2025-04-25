"""
Main topology generation module.

This module handles the generation of molecular topology files using different force fields
(AMBER99 and MARTINI3). It provides functionality for file management, cleanup, and 
topology file organization.
"""

import os
from pathlib import Path
from typing import Any, Optional, Set, List
import shutil
from colorama import init, Fore, Style

from colbuilder.core.utils.files import FileManager, managed_resources
from colbuilder.core.geometry.system import System
from colbuilder.core.topology.amber import Amber, build_amber99
from colbuilder.core.topology.martini import Martini, build_martini3
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import (
    TopologyGenerationError
)
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

REQUIRED_FF_FILES: List[str] = ['residuetypes.dat', 'specbond.dat']
TEMP_FILES_TO_CLEAN: Set[str] = {
    'col_*.*_go-*.itp', 
    'col_*.*.itp', 
    'col_[0-9]*_go-sites.itp',
    'col_[0-9]*_go-table.itp',
    'col_[0-9]*_go-harm.itp',
    'go-table.itp',
    'col_go-sites.itp',
    '*.CG.pdb',
    'D'
}

def cleanup_temporary_files(ff_name: str, temp_patterns: Set[str], search_dirs: Optional[List[Path]] = None) -> None:
    """
    Remove temporary files and directories created during topology generation.
    
    Parameters
    ----------
    ff_name : str
        Force field directory name to be cleaned up
    temp_patterns : Set[str]
        Patterns matching temporary files and directories to remove
    search_dirs : Optional[List[Path]]
        List of specific directories to search in. If None, searches in current directory
    """
    try:
        # Use current directory if no search dirs specified
        dirs_to_search = search_dirs if search_dirs else [Path()]
        
        for search_dir in dirs_to_search:
            if not search_dir.exists():
                LOG.warning(f"Search directory does not exist: {search_dir}")
                continue
                
            for pattern in temp_patterns:
                if "*" not in pattern:
                    path_obj = search_dir / pattern
                    if path_obj.exists():
                        if path_obj.is_dir():
                            shutil.rmtree(path_obj)
                        else:
                            os.remove(path_obj)
                else:
                    for matched_path in search_dir.glob(pattern):
                        if matched_path.is_dir():
                            shutil.rmtree(matched_path)
                        else:
                            os.remove(matched_path)
        
        # Handle force field directory cleanup
        for search_dir in dirs_to_search:
            copied_ff_dir = search_dir / ff_name
            if copied_ff_dir.exists():
                shutil.rmtree(copied_ff_dir)
            
    except Exception as e:
        LOG.warning(f"Error during cleanup: {str(e)}")


def setup_topology_directory(system_name: str, ff_name: str) -> Path:
    """
    Create and prepare a directory for topology file generation.
    
    Parameters
    ----------
    system_name : str
        System identifier used for directory naming
        
    Returns
    -------
    Path
        Path to the created topology directory
    """
    topology_dir = Path(f"{system_name}_{ff_name}_topology_files")
    topology_dir.mkdir(exist_ok=True)
    return topology_dir


def organize_topology_files(topology_dir: Path, species: str) -> None:
    """
    Organize and move topology files to their final location.
    
    Parameters
    ----------
    topology_dir : Path
        Destination directory for topology files
    species : str
        Species identifier used in file naming
    """
    try:
        for top_file in Path().glob(f"collagen_fibril_*.top"):
            shutil.copy2(top_file, topology_dir / top_file.name)

        for gro_file in Path().glob("*.gro"):
            shutil.copy2(gro_file, topology_dir / gro_file.name)
            
        for itp_file in Path().glob("*.itp"):
            shutil.copy2(itp_file, topology_dir / itp_file.name)
            os.remove(itp_file)  
            
    except Exception as e:
        LOG.warning(f"Error organizing topology files: {str(e)}")


@timeit
async def build_topology(system: System, config: ColbuilderConfig, file_manager: Optional[FileManager] = None) -> Any:
    """
    Build molecular system topology and organize output files.
    
    Handles the complete topology generation process including:
    - Setting up working directories
    - Copying necessary geometry files
    - Generating topology based on selected force field
    - Organizing output files
    
    Parameters
    ----------
    system : System
        Molecular system to process
    config : ColbuilderConfig
        Configuration settings including force field parameters
    file_manager : Optional[FileManager]
        File manager for handling I/O operations
        
    Returns
    -------
    Any
        Processed molecular system
        
    Raises
    ------
    TopologyGenerationError
        If topology generation fails at any stage
    """
    try:
        file_manager = file_manager or FileManager(config)
        topology_dir = file_manager.get_temp_dir("topology_gen")
        original_dir = Path.cwd()
        
        try:
            os.chdir(topology_dir)
            geometry_dir = Path(".tmp/geometry_gen")
            if not geometry_dir.exists():
                geometry_dir = Path(original_dir) / ".tmp" / "geometry_gen"
            
            if geometry_dir.exists():
                cap_files = list(geometry_dir.glob("**/*.caps.pdb"))
                
                if list(system.get_models()):
                    first_model = system.get_model(model_id=list(system.get_models())[0])
                    model_type = first_model.type
                    type_dir = topology_dir / model_type
                    type_dir.mkdir(exist_ok=True, parents=True)
                    
                    for cap_file in cap_files:
                        dest_file = type_dir / cap_file.name
                        shutil.copy(cap_file, dest_file)
            
            force_field = config.force_field
            
            if force_field == 'amber99':
                ff = f"{force_field}sb-star-ildnp"
                LOG.subsection(f'Building topology based on the {force_field} force field')
                await build_amber99(system=system, config=config, file_manager=file_manager)
                
            elif force_field == 'martini3':
                ff = force_field
                LOG.subsection(f'Building topology based on the {force_field} force field')
                
                from colbuilder.core.utils.martinize_finder import find_and_install_custom_force_field
                find_and_install_custom_force_field(config.FORCE_FIELD_DIR)
                
                await build_martini3(system=system, config=config, file_manager=file_manager)
            else:
                raise TopologyGenerationError(
                    message='Invalid or unsupported force field specified',
                    error_code="TOP_ERR_001",
                    context={"force_field": force_field}
                )
                
            output_topology_dir = file_manager.ensure_dir(f"{config.species}_topology_files")
            
            for top_file in Path().glob(f"collagen_fibril_*.top"):
                file_manager.copy_to_directory(top_file, dest_dir=output_topology_dir)
                
            for itp_file in Path().glob("col_[0-9]*.itp"):
                file_manager.copy_to_directory(itp_file, dest_dir=output_topology_dir)
                    
            LOG.info(f"{Fore.BLUE}Topology files written at: {output_topology_dir}{Style.RESET_ALL}")
                
            return system
            
        finally:
            os.chdir(original_dir)
            search_dirs = [topology_dir, output_topology_dir]
            if not config.topology_debug:
                cleanup_temporary_files(config.force_field, TEMP_FILES_TO_CLEAN, search_dirs=search_dirs)
            
    except TopologyGenerationError:
        raise
    except Exception as e:
        raise TopologyGenerationError(
            message="Unexpected error in topology generation",
            original_error=e,
            error_code="TOP_ERR_001",
            context={"force_field": config.force_field}
        )