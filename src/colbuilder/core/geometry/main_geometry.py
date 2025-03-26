"""
Colbuilder Main Geometry Module

Main entry point for geometry operations, coordinating the crystal,
mixing, and replacement services.
"""

import os
import time
import shutil
import traceback
from pathlib import Path
from typing import Optional, Set, List, Union

from ..utils.exceptions import GeometryGenerationError
from ..utils.config import ColbuilderConfig
from ..utils.logger import setup_logger
from .crystal_builder import CrystalBuilder
from .crosslink_mixer import CrosslinkMixer
from .geometry_replacer import CrosslinkReplacer
from .system import System
from .crystal import Crystal
from colorama import init, Fore, Style

LOG = setup_logger(__name__)

# Standard temporary files and directories created during processing
STANDARD_TEMP_FILES = {"replace.txt"}
STANDARD_TEMP_DIRS = {}
# STANDARD_TEMP_DIRS = {"NC", "NCP", "D", "T"}

def cleanup_temp_files(
    temp_files: Optional[Set[str]] = None, 
    temp_dirs: Optional[Set[str]] = None,
    include_standard: bool = True
) -> None:
    """
    Clean up temporary files and directories.
    
    Args:
        temp_files: Optional set of specific temporary files to clean up
        temp_dirs: Optional set of specific temporary directories to clean up
        include_standard: Whether to also clean up standard temporary files/directories
    """
    files_to_clean = set(temp_files) if temp_files else set()
    dirs_to_clean = set(temp_dirs) if temp_dirs else set()
    
    if include_standard:
        files_to_clean.update(STANDARD_TEMP_FILES)
        dirs_to_clean.update(STANDARD_TEMP_DIRS)
    
    # Clean up directories
    for dir_path in dirs_to_clean:
        if os.path.exists(dir_path) and os.path.isdir(dir_path):
            try:
                shutil.rmtree(dir_path)
                LOG.debug(f"Removed temporary directory: {dir_path}")
            except Exception as e:
                LOG.warning(f"Failed to remove temporary directory {dir_path}: {str(e)}")
    
    # Clean up files
    for file_path in files_to_clean:
        if os.path.exists(file_path) and os.path.isfile(file_path):
            try:
                os.remove(file_path)
                LOG.debug(f"Removed temporary file: {file_path}")
            except Exception as e:
                LOG.warning(f"Failed to remove temporary file {file_path}: {str(e)}")

class GeometryService:
    """
    Main service coordinating geometry operations.
    
    This class orchestrates the different components involved in collagen
    geometry generation, mixing, and crosslink replacement.
    """
    
    def __init__(self, config: ColbuilderConfig) -> None:
        """
        Initialize geometry service with configuration.
        
        Args:
            config: Configuration settings for geometry operations
        """
        self.config = config
        self.crystal_service = CrystalBuilder()
        self.mixer_service = CrosslinkMixer() 
        self.replacer_service = CrosslinkReplacer()  
        self.temp_files: Set[str] = set()
        self.temp_dirs: Set[str] = set()
        
    def _track_temp_resources(self, files: Optional[List[str]] = None, dirs: Optional[List[str]] = None) -> None:
        """
        Track temporary files and directories for later cleanup.
        
        Args:
            files: List of file paths to track
            dirs: List of directory paths to track
        """
        if files:
            self.temp_files.update(files)
            for file_path in files:
                LOG.debug(f"Tracking temporary file: {file_path}")
                
        if dirs:
            self.temp_dirs.update(dirs)
            for dir_path in dirs:
                LOG.debug(f"Tracking temporary directory: {dir_path}")
    
    def _cleanup(self) -> None:
        """Clean up all tracked temporary resources."""
        LOG.info("Cleaning up temporary files...")
        cleanup_temp_files(self.temp_files, self.temp_dirs)
        
    async def _handle_direct_replacement(self) -> Optional[System]:
        """
        Handle direct replacement without geometry generation.
        
        Returns:
            None as direct replacement doesn't produce a System object
            
        Raises:
            GeometryGenerationError: If replacement fails
        """
        LOG.info("Starting direct replacement mode")
        
        # Track standard temp resources
        self._track_temp_resources(
            files=list(STANDARD_TEMP_FILES),
            dirs=list(STANDARD_TEMP_DIRS)
        )
        
        try:
            await self.replacer_service.replace_direct(self.config)
            LOG.info("Direct replacement completed successfully")
            return None
        except Exception as e:
            LOG.error(f"Error during direct replacement: {str(e)}")
            raise GeometryGenerationError(
                message=f"Failed to complete direct replacement: {str(e)}",
                original_error=e,
                error_code="GEO_ERR_004"
            )
    
    async def _handle_mixing_only(self) -> System:
        """
        Handle mixing operation without geometry generation.
        
        Returns:
            Mixed system
            
        Raises:
            GeometryGenerationError: If mixing fails
        """
        if not self.config.files_mix:
            raise GeometryGenerationError(
                message="No files_mix provided for mixing operation",
                error_code="GEO_ERR_003"
            )
        
        # Create a basic system from the first file in files_mix
        system = System()
        crystal = Crystal(pdb=str(self.config.files_mix[0]))
        system.crystal = crystal
        
        # Apply mixing
        system = await self.mixer_service.mix(system, self.config)
        
        # Write the final mixed system to PDB
        output_pdb = f"{self.config.output}.pdb"
        LOG.info(f"Writing mixed system to {output_pdb}")
        system.write_pdb(pdb_out=self.config.output, fibril_length=self.config.fibril_length)
        
        return system
    
    async def _handle_full_generation(self) -> Optional[System]:
        """
        Handle full geometry generation with optional mixing and replacement.
        
        Returns:
            Generated system
            
        Raises:
            GeometryGenerationError: If any step fails
        """
        system = None
        
        # Generate base geometry if requested
        if self.config.geometry_generator:
            LOG.info('Geometry mode: building system from crystal')
            system = await self.crystal_service.build(self.config)
        
        # Apply mixing if requested
        if self.config.mix_bool and system:
            LOG.section("Mixing geometry...")
            LOG.subsection("Mixing Geometry")
            system = await self.mixer_service.mix(system, self.config)
            LOG.info(f"{Fore.BLUE}Mixing completed.{Style.RESET_ALL}")
        
        # Apply replacement if requested
        if self.config.replace_bool and system:
            LOG.section("Replacing crosslinks...")
            LOG.subsection("Replacing Geometry")
            
            # Track standard temp resources
            self._track_temp_resources(
                files=list(STANDARD_TEMP_FILES),
                dirs=list(STANDARD_TEMP_DIRS)
            )
            
            # Allow time for file operations to complete
            time.sleep(1)
            
            # Apply the replacement
            system = await self.replacer_service.replace_in_system(system, self.config)
            LOG.info(f"{Fore.BLUE}Replacement completed.{Style.RESET_ALL}")
        
        # Write the final output PDB if we have a system
        if system:
            output_pdb = f"{self.config.output}.pdb"
            LOG.info(f"Writing final system to {output_pdb}")
            
            # Don't let System.write_pdb clean up when replacement is enabled
            # to avoid conflicts with our own cleanup
            # cleanup = not self.config.replace_bool
            system.write_pdb(
                pdb_out=self.config.output, 
                fibril_length=self.config.fibril_length,
                cleanup=False
            )
        
        return system
        
    async def build_geometry(self) -> Optional[System]:
        """
        Main entry point for geometry generation.
        
        This method orchestrates the entire geometry generation process based on
        the configuration settings, handling crystal building, mixing, and
        crosslink replacement.
        
        Returns:
            Generated system or None if only direct replacement was performed
            
        Raises:
            GeometryGenerationError: If any step in the process fails
        """
        try:
            # Handle direct replacement case (no geometry generation)
            if self.config.replace_bool and not self.config.geometry_generator:
                return await self._handle_direct_replacement()
            
            # Handle mixing-only case (no geometry generation)
            if self.config.mix_bool and not self.config.geometry_generator:
                return await self._handle_mixing_only()
            
            # Handle informational case (no operations requested)
            if not self.config.geometry_generator and not self.config.replace_bool and not self.config.mix_bool:
                LOG.info('Set -geometry flag to generate microfibrillar structure PDB file')
                return None
            
            # Handle full generation case (with optional mixing and replacement)
            return await self._handle_full_generation()
            
        except GeometryGenerationError:
            # Just re-raise GeometryGenerationError as it's already properly formatted
            raise
        except Exception as e:
            # Convert any other exception to GeometryGenerationError
            LOG.error(f"Unexpected error in geometry generation: {str(e)}")
            traceback.print_exc()
            raise GeometryGenerationError(
                message=f"Unexpected error in geometry generation: {str(e)}",
                original_error=e,
                error_code="GEO_ERR_001", 
                context={"config": self.config.model_dump()}
            )
        # finally:
        #     if not self.config.topology_generator:
        #         self._cleanup()

async def build_geometry(config: ColbuilderConfig) -> Optional[System]:
    """
    Build geometry from configuration.
    
    Standalone function that creates a GeometryService and invokes its
    build_geometry method.
    
    Args:
        config: Configuration for geometry operations
        
    Returns:
        Generated system or None
    """
    service = GeometryService(config)
    return await service.build_geometry()

async def mix_geometry(system: System, config: ColbuilderConfig) -> System:
    """
    Mix geometry types in system.
    
    Standalone function that creates a CrosslinkMixer and invokes its
    mix method.
    
    Args:
        system: System containing models to be mixed
        config: Configuration for mixing operation
        
    Returns:
        Mixed system
    """
    mixer = CrosslinkMixer()
    return await mixer.mix(system, config)

async def replace_geometry(system: System, config: ColbuilderConfig) -> Optional[System]:
    """
    Replace crosslinks in system according to replacement ratio.
    
    Standalone function that handles both direct replacement (when no system
    is provided) and system-based replacement.
    
    Args:
        system: System containing models with crosslinks, or None for direct replacement
        config: Configuration for replacement
        
    Returns:
        System with replaced crosslinks, or None for direct replacement
    """
    # Track temporary resources for cleanup
    temp_files = set(STANDARD_TEMP_FILES)
    temp_dirs = set(STANDARD_TEMP_DIRS)
    
    try:
        replacer = CrosslinkReplacer()
        
        # Determine if this is a direct replacement case
        is_direct_replacement = (
            system is None and 
            config.replace_file and 
            os.path.exists(config.replace_file) and
            _is_pdb_file(config.replace_file)
        )
        
        if is_direct_replacement:
            # Direct replacement path
            LOG.info("Using direct replacement approach")
            await replacer.replace_direct(config)
            LOG.info(f"PDB file written to {config.output}.pdb")
            
            # Try to create a minimal system for consistency
            output_pdb = f"{config.output}.pdb"
            if os.path.exists(output_pdb):
                try:
                    crystal = Crystal(pdb=output_pdb)
                    return System(crystal=crystal)
                except Exception as e:
                    LOG.warning(f"Could not create minimal system after direct replacement: {e}")
            
            return None
        else:
            # System-based replacement path
            if not system:
                raise GeometryGenerationError(
                    message="No system provided for replacement and input file is not a valid PDB",
                    error_code="GEO_ERR_004"
                )
                
            system = await replacer.replace_in_system(system, config)
            
            # Write the output PDB file
            output_pdb = f"{config.output}.pdb"
            system.write_pdb(pdb_out=config.output, fibril_length=config.fibril_length, cleanup=False)
            LOG.info(f"PDB file written to {output_pdb}")
            
            return system
    finally:
        cleanup_temp_files(temp_files)

def _is_pdb_file(file_path: str) -> bool:
    """
    Check if a file is a PDB file based on its contents.
    
    Args:
        file_path: Path to the file to check
        
    Returns:
        True if the file appears to be a PDB file, False otherwise
    """
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            return first_line.startswith(("ATOM", "CRYST1", "HETATM"))
    except Exception:
        return False