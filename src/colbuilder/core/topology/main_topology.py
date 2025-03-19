# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import subprocess
import os
from pathlib import Path
from typing import Any, Optional, Set
import shutil
import asyncio
from colorama import init, Fore, Style

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
TEMP_FILES_TO_CLEAN = ['*.itp', '*.CG.pdb', '*.merge.pdb', 'create_goVirt.py', 'tmp.pdb', '*.top', 'map.*', 'contactmap', 'amber99*', '*.dat', 'D', 'T', 'NC']


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
async def build_topology(system: System, config: ColbuilderConfig) -> Any:
    """
    Build the topology of a molecular system and organize output files.
    
    Parameters
    ----------
    system : System
        The molecular system to process
    config : ColbuilderConfig
        Configuration object containing settings
        
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
        force_field = config.force_field
        
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            LOG.debug(f"Model {model_id} - Type: {model.type}, Connect: {model.connect}")
            
        if force_field == 'amber99':
            ff = f"{force_field}sb-star-ildnp"
            LOG.subsection(f'Building topology based on the {force_field} force field')
            await build_amber99(system=system, config=config)
            
        elif force_field == 'martini3':
            ff = f"{force_field}"
            LOG.subsection(f'Building topology based on the {force_field} force field')
            
            # Install custom force field files before building the topology
            LOG.info("Installing Martini 3.0 force field custom files")
            from colbuilder.core.utils.martinize_finder import find_and_install_custom_force_field
            
            if not find_and_install_custom_force_field(config.FORCE_FIELD_DIR):
                LOG.warning("Failed to install custom force field files, proceeding with existing installation")
            
            await build_martini3(system=system, config=config)
        else:
            raise TopologyGenerationError(
                message='Invalid or unsupported force field specified',
                error_code="TOP_ERR_001",
                context={"force_field": force_field}
            )
            
        # Organize topology files in the final top directory
        LOG.debug("Setting up topology directory and organizing files")
        topology_dir = setup_topology_directory(config.species)
        organize_topology_files(topology_dir, config.species)
        
        # Clean up temporary files
        LOG.debug("Cleaning up temporary files")
        cleanup_temporary_files(force_field, TEMP_FILES_TO_CLEAN)
            
        return system
    except TopologyGenerationError:
        try:
            LOG.warning("Attempting cleanup after error")
            cleanup_temporary_files(config.force_field, TEMP_FILES_TO_CLEAN)
        except Exception as cleanup_error:
            LOG.warning(f"Cleanup after error failed: {cleanup_error}")
        raise
    except Exception as e:
        try:
            LOG.warning("Attempting cleanup after unexpected error")
            cleanup_temporary_files(config.force_field, TEMP_FILES_TO_CLEAN)
        except Exception as cleanup_error:
            LOG.warning(f"Cleanup after error failed: {cleanup_error}")
        
        raise TopologyGenerationError(
            message="Unexpected error in topology generation",
            original_error=e,
            error_code="TOP_ERR_001",
            context={"force_field": config.force_field}
        )