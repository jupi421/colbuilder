# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import subprocess
import os
from pathlib import Path
from typing import Any, Optional, Set
import shutil
import traceback
import asyncio
from colorama import init, Fore, Style

from ..utils.files import FileManager, managed_resources
from colbuilder.core.geometry.system import System
from colbuilder.core.topology.amber import Amber, build_amber99
from colbuilder.core.topology.martini import Martini, build_martini3
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import (
    TopologyGenerationError,
    ColbuilderError,
    ErrorCategory,
    ErrorSeverity,
    ColbuilderErrorDetail
)
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

REQUIRED_FF_FILES = ['residuetypes.dat', 'specbond.dat']
TEMP_FILES_TO_CLEAN = ['*.itp', '*.CG.pdb', '*.merge.pdb', 'create_goVirt.py', 'tmp.pdb', '*.top', 'map.*', 'contactmap', 'amber99*', '*.dat']


def cleanup_temporary_files(ff_name: str, temp_patterns: Set[str]) -> None:
    """
    Clean up temporary files and directories.
    
    Parameters
    ----------
    ff_name : str
        Name of the force field directory
    temp_patterns : Set[str]
        Set of patterns for temporary files and directories to remove
    """
    try:
        for pattern in temp_patterns:
            if not "*" in pattern:
                path_obj = Path(pattern)
                if path_obj.exists():
                    if path_obj.is_dir():
                        shutil.rmtree(path_obj)
                        LOG.debug(f"Removed temporary directory: {path_obj}")
                    else:
                        os.remove(path_obj)
                        LOG.debug(f"Removed temporary file: {path_obj}")
            else:
                for matched_path in Path().glob(pattern):
                    if matched_path.is_dir():
                        shutil.rmtree(matched_path)
                        LOG.debug(f"Removed temporary directory: {matched_path}")
                    else:
                        os.remove(matched_path)
                        LOG.debug(f"Removed temporary file: {matched_path}")
                
        copied_ff_dir = Path(ff_name)
        if copied_ff_dir.exists():
            shutil.rmtree(copied_ff_dir)
            LOG.debug(f"Removed copied force field directory: {copied_ff_dir}")
    except Exception as e:
        LOG.warning(f"Error during cleanup: {str(e)}")


def setup_topology_directory(system_name: str) -> Path:
    """
    Create and return path to topology directory.
    
    Parameters
    ----------
    system_name : str
        Name of the system for the directory
        
    Returns
    -------
    Path
        Path to the created topology directory
    """
    topology_dir = Path(f"{system_name}_topology_files")
    topology_dir.mkdir(exist_ok=True)
    return topology_dir


def organize_topology_files(topology_dir: Path, species: str) -> None:
    """
    Move topology files to the final directory.
    
    Parameters
    ----------
    topology_dir : Path
        Path to topology directory
    species : str
        Species name for naming convention
    """
    try:
        # Copy topology files
        for top_file in Path().glob(f"collagen_fibril_*.top"):
            shutil.copy2(top_file, topology_dir / top_file.name)
            
        # Copy and remove ITP files
        for itp_file in Path().glob("*.itp"):
            shutil.copy2(itp_file, topology_dir / itp_file.name)
            os.remove(itp_file)  
            
        LOG.info(f"{Fore.BLUE}Topology files written at: {topology_dir}{Style.RESET_ALL}")
    except Exception as e:
        LOG.warning(f"Error organizing topology files: {str(e)}")


@timeit
async def build_topology(system: System, config: ColbuilderConfig, file_manager: Optional[FileManager] = None) -> Any:
    """
    Build the topology of a molecular system and organize output files.
    
    Parameters
    ----------
    system : System
        The molecular system to process
    config : ColbuilderConfig
        Configuration object containing settings
    file_manager : Optional[FileManager]
        File manager for consistent file handling
        
    Returns
    -------
    Any
        The processed system
        
    Raises
    ------
    TopologyGenerationError
        If topology generation fails
    """
    try:
        # Initialize file manager if not provided
        if file_manager is None:
            file_manager = FileManager(config)
        
        # Get the topology directory from the file manager
        topology_dir = file_manager.get_temp_dir("topology_gen")
        original_dir = Path.cwd()
        
        try:
            # Change to topology directory
            os.chdir(topology_dir)
            
            # Search for cap files and copy them to the topology directory
            geometry_dir = Path(".tmp/geometry_gen")
            if not geometry_dir.exists():
                geometry_dir = Path(original_dir) / ".tmp" / "geometry_gen"
            
            if geometry_dir.exists():
                cap_files = list(geometry_dir.glob("**/*.caps.pdb"))
                
                # Get the model type from the first model
                if len(list(system.get_models())) > 0:
                    first_model = system.get_model(model_id=list(system.get_models())[0])
                    model_type = first_model.type
                    type_dir = topology_dir / model_type
                    type_dir.mkdir(exist_ok=True, parents=True)
                    
                    # Copy cap files to the topology directory
                    for cap_file in cap_files:
                        dest_file = type_dir / cap_file.name
                        shutil.copy(cap_file, dest_file)
                else:
                    LOG.warning("No models found in system, cannot determine model type")
            else:
                LOG.warning(f"Geometry directory not found: {geometry_dir}")
            
            force_field = config.force_field
            
            for model_id in system.get_models():
                model = system.get_model(model_id=model_id)
                LOG.debug(f"Model {model_id} - Type: {model.type}, Connect: {model.connect}")
                
            if force_field == 'amber99':
                ff = f"{force_field}sb-star-ildnp"
                LOG.subsection(f'Building topology based on the {force_field} force field')
                await build_amber99(system=system, config=config, file_manager=file_manager)
                
            elif force_field == 'martini3':
                ff = f"{force_field}"
                LOG.subsection(f'Building topology based on the {force_field} force field')
                
                # Install custom force field files before building the topology
                LOG.info("Installing Martini 3.0 force field custom files")
                from colbuilder.core.utils.martinize_finder import find_and_install_custom_force_field
                
                if not find_and_install_custom_force_field(config.FORCE_FIELD_DIR):
                    LOG.warning("Failed to install custom force field files, proceeding with existing installation")
                
                await build_martini3(system=system, config=config, file_manager=file_manager)
            else:
                raise TopologyGenerationError(
                    message='Invalid or unsupported force field specified',
                    error_code="TOP_ERR_001",
                    context={"force_field": force_field}
                )
                
            # Organize topology files in the final top directory
            LOG.debug("Setting up topology directory and organizing files")
            output_topology_dir = file_manager.ensure_dir(f"{config.species}_topology_files")
            
            # Copy topology files
            for top_file in Path().glob(f"collagen_fibril_*.top"):
                dest_path = file_manager.copy_to_directory(top_file, dest_dir=output_topology_dir)
                LOG.debug(f"Copied topology file to: {dest_path}")
                
            # Copy ITP files
            for itp_file in Path().glob("*.itp"):
                dest_path = file_manager.copy_to_directory(itp_file, dest_dir=output_topology_dir)
                LOG.debug(f"Copied ITP file to: {dest_path}")
                    
            LOG.info(f"{Fore.BLUE}Topology files written at: {output_topology_dir}{Style.RESET_ALL}")
                
            return system
            
        finally:
            # Return to original directory
            os.chdir(original_dir)
            LOG.debug(f"Returned to original directory: {original_dir}")
            
            # Clean up temporary files (unless in debug mode)
            if not config.debug:
                cleanup_temporary_files(config.force_field, TEMP_FILES_TO_CLEAN)
            
    except TopologyGenerationError:
        raise
    except Exception as e:
        LOG.error(f"Unexpected error in topology generation: {str(e)}")
        LOG.debug(f"Traceback: {traceback.format_exc()}")
        
        raise TopologyGenerationError(
            message="Unexpected error in topology generation",
            original_error=e,
            error_code="TOP_ERR_001",
            context={"force_field": config.force_field}
        )