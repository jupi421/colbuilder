"""
Colbuilder Main Geometry Module

Main entry point for geometry operations, coordinating the crystal,
mixing, and replacement services.
"""

import os
import time
import shutil
import traceback
import asyncio
from pathlib import Path
from typing import Optional, Set, List, Union, Tuple

from ..utils.exceptions import GeometryGenerationError
from ..utils.config import ColbuilderConfig
from ..utils.logger import setup_logger
from ..utils.files import FileManager, managed_resources
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
                LOG.info(f"Removed temporary directory: {dir_path}")
            except Exception as e:
                LOG.warning(f"Failed to remove temporary directory {dir_path}: {str(e)}")
    
    # Clean up files
    for file_path in files_to_clean:
        if os.path.exists(file_path) and os.path.isfile(file_path):
            try:
                os.remove(file_path)
                LOG.info(f"Removed temporary file: {file_path}")
            except Exception as e:
                LOG.warning(f"Failed to remove temporary file {file_path}: {str(e)}")

class GeometryService:
    """
    Main service coordinating geometry operations.
    
    This class orchestrates the different components involved in collagen
    geometry generation, mixing, and crosslink replacement.
    """
    
    def __init__(self, config: ColbuilderConfig, file_manager: Optional[FileManager] = None) -> None:
        """
        Initialize geometry service with configuration.
        
        Args:
            config: Configuration settings for geometry operations
            file_manager: Optional file manager for consistent file handling
        """
        self.config = config
        self.crystal_service = CrystalBuilder()
        self.mixer_service = CrosslinkMixer() 
        self.replacer_service = CrosslinkReplacer()  
        self.temp_files: Set[str] = set()
        self.temp_dirs: Set[str] = set()
        self.file_manager = file_manager or FileManager(config)
        self.original_dir = Path.cwd()
        self.temp_dir = None
        
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
                LOG.info(f"Tracking temporary file: {file_path}")
                
        if dirs:
            self.temp_dirs.update(dirs)
            for dir_path in dirs:
                LOG.info(f"Tracking temporary directory: {dir_path}")
    
    def _cleanup(self) -> None:
        """Clean up all tracked temporary resources."""
        LOG.info("Cleaning up temporary files...")
        cleanup_temp_files(self.temp_files, self.temp_dirs)
        
    @managed_resources("geometry_operation")
    async def _handle_mixing_only(self) -> Tuple[Path, Path]:
        """
        Handle mixing operation without geometry generation.
        
        Returns:
            Tuple containing temporary directory path and output PDB file path
                
        Raises:
            GeometryGenerationError: If mixing fails
        """
        import traceback
        
        if not self.config.files_mix:
            raise GeometryGenerationError(
                message="No files_mix provided for mixing operation",
                error_code="GEO_ERR_003"
            )
        
        # Create a temporary directory for mixing operations
        temp_dir = self.file_manager.get_temp_path("mix_operation", create_dir=True)
        LOG.info(f"{Fore.BLUE}Created temporary directory for mixing: {temp_dir}{Style.RESET_ALL}")

        # Store original directory to restore later
        original_dir = Path.cwd()
        
        try:
            # Parse ratio_mix if it's a string
            if isinstance(self.config.ratio_mix, str):
                LOG.info(f"Converting ratio_mix string to dictionary: {self.config.ratio_mix}")
                ratio_dict = {}
                for part in self.config.ratio_mix.split():
                    if ':' in part:
                        key, value = part.split(':')
                        try:
                            ratio_dict[key] = int(value)
                        except ValueError:
                            LOG.error(f"Invalid ratio value in {part}")
                            ratio_dict[key] = 0
                self.config.ratio_mix = ratio_dict
                LOG.info(f"Converted ratio_mix: {self.config.ratio_mix}")
            
            # Ensure ratio_mix is a dictionary
            if not isinstance(self.config.ratio_mix, dict):
                LOG.error(f"ratio_mix is not a dictionary: {type(self.config.ratio_mix)}")
                raise GeometryGenerationError(
                    message=f"Invalid ratio_mix format: {self.config.ratio_mix}",
                    error_code="GEO_ERR_009"
                )
            
            # Copy mix files to the temporary directory
            updated_mix_files = []
            for mix_file in self.config.files_mix:
                try:
                    LOG.debug(f"Copying file {mix_file} to {temp_dir}")
                    dest_path = self.file_manager.copy_to_directory(mix_file, dest_dir=temp_dir)
                    updated_mix_files.append(dest_path)
                    LOG.debug(f"File copied successfully to {dest_path}")
                except Exception as e:
                    LOG.error(f"Failed to copy file {mix_file}: {str(e)}")
                    raise GeometryGenerationError(
                        message=f"Failed to copy file {mix_file}: {str(e)}",
                        error_code="GEO_ERR_010"
                    )
            
            # Check that we have the right number of files
            if len(updated_mix_files) < len(self.config.ratio_mix):
                LOG.warning(f"Not enough files ({len(updated_mix_files)}) for ratio_mix types ({len(self.config.ratio_mix)})")
            
            # Update config with local paths
            original_files_mix = self.config.files_mix
            self.config.files_mix = updated_mix_files
            
            # Create a basic system from the first file in files_mix
            LOG.debug(f"Creating system from first file: {self.config.files_mix[0]}")
            system = System()
            crystal = Crystal(pdb=str(self.config.files_mix[0]))
            system.crystal = crystal
            
            # Apply mixing with temporary directory
            LOG.debug(f"Calling mixer_service.mix with temp_dir: {temp_dir}")
            
            # Ensure temp_dir is properly set
            if not isinstance(temp_dir, Path):
                LOG.error(f"temp_dir is not a Path object: {type(temp_dir)}")
                raise GeometryGenerationError(
                    message="Invalid temporary directory type",
                    error_code="GEO_ERR_007"
                )
                
            if not temp_dir.exists():
                LOG.error(f"temp_dir does not exist: {temp_dir}")
                raise GeometryGenerationError(
                    message="Temporary directory does not exist",
                    error_code="GEO_ERR_008"
                )
            
            # Call the mixer service
            try:
                system = await self.mixer_service.mix(system, self.config, temp_dir)
            except Exception as e:
                LOG.error(f"Mixer service failed: {str(e)}")
                LOG.debug(f"Traceback: {traceback.format_exc()}")
                raise GeometryGenerationError(
                    message=f"Mixer service failed: {str(e)}",
                    error_code="GEO_ERR_011",
                    original_error=e
                )
            
            # Write the final mixed system to PDB in the output location
            try:
                output_pdb = self.file_manager.get_output_path(self.config.output, ".pdb")
                LOG.info(f"{Fore.BLUE}Writing mixed system to {output_pdb}{Style.RESET_ALL}")
                
                # Pass the temporary directory to write_pdb to help find caps files
                system.write_pdb(pdb_out=output_pdb, fibril_length=self.config.fibril_length, temp_dir=temp_dir)
            except Exception as e:
                LOG.error(f"Failed to write output PDB: {str(e)}")
                raise GeometryGenerationError(
                    message=f"Failed to write output PDB: {str(e)}",
                    error_code="GEO_ERR_012",
                    original_error=e
                )
            
            # Return both temp_dir and output_pdb paths for consistency with other handlers
            return temp_dir, output_pdb
            
        except GeometryGenerationError:
            # Re-raise GeometryGenerationError directly
            raise
        except Exception as e:
            error_msg = f"Mixing operation failed: {str(e)}"
            LOG.error(error_msg)
            LOG.debug(f"Traceback: {traceback.format_exc()}")
            raise GeometryGenerationError(
                message=error_msg,
                error_code="GEO_ERR_005",
                original_error=e
            )
            
        finally:
            # Restore original files_mix in config
            if 'original_files_mix' in locals():
                self.config.files_mix = original_files_mix
                
            # Return to original directory
            os.chdir(original_dir)

    @managed_resources("geometry_operation")
    async def _handle_full_generation(self) -> Tuple[Optional[System], Optional[Path]]:
        """
        Handle full geometry generation with optional mixing and replacement.

        Returns:
            Generated system

        Raises:
            GeometryGenerationError: If any step fails
        """
        try:
            # Create the temporary directory for geometry generation
            if not self.temp_dir:
                self.temp_dir = self.file_manager.get_temp_path("geometry_gen", create_dir=True)
            LOG.info(f"Using temporary directory for geometry generation: {self.temp_dir}")

            # Change to the temporary directory
            os.chdir(self.temp_dir)

            temp_config = self.config.copy()
            if self.config.pdb_file:
                pdb_path = Path(self.config.pdb_file)
                if pdb_path.exists():
                    local_pdb = self.temp_dir / pdb_path.name
                    shutil.copy2(pdb_path, local_pdb)
                    temp_config.pdb_file = local_pdb
                else:
                    raise GeometryGenerationError(
                        message=f"PDB file not found: {pdb_path}",
                        error_code="GEO_ERR_001"
                    )

            temp_config.working_directory = self.temp_dir

            system = None

            # Perform geometry generation
            if temp_config.geometry_generator:
                LOG.info('Geometry mode: building system from crystal')
                system = await self.crystal_service.build(temp_config)

            # Perform mixing if required
            if temp_config.mix_bool and system:
                LOG.section("Mixing geometry...")
                system = await self.mixer_service.mix(system, temp_config)
                LOG.info("Mixing completed.")

            # Perform replacement if required
            if temp_config.replace_bool and system:
                LOG.section("Replacing crosslinks...")
                system = await self.replacer_service.replace_in_system(system, temp_config)
                LOG.info("Replacement completed.")

            # Write the final system to the output file
            if system:
                output_prefix = temp_config.species
                if temp_config.output:
                    output_prefix = temp_config.output

                output_pdb_path = self.temp_dir / f"{output_prefix}.pdb"
                LOG.info(f"Writing final system to {output_pdb_path}")

                system.write_pdb(
                    pdb_out=output_prefix,
                    fibril_length=temp_config.fibril_length,
                    cleanup=False
                )

                final_pdb = self.file_manager.copy_to_output(output_pdb_path)

                return system, final_pdb

            return None, None
        except Exception as e:
            LOG.error(f"Error during full generation: {str(e)}")
            raise GeometryGenerationError(
                message=f"Failed to complete geometry generation: {str(e)}",
                original_error=e,
                error_code="GEO_ERR_001"
            )
        finally:
            os.chdir(self.original_dir)
        
    async def build_geometry(self) -> Tuple[Optional[System], Optional[Path]]:
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
            if self.config.replace_bool and not self.config.geometry_generator:
                return await self._handle_direct_replacement()
            
            if self.config.mix_bool and not self.config.geometry_generator:
                return await self._handle_mixing_only()
            
            if not self.config.geometry_generator and not self.config.replace_bool and not self.config.mix_bool:
                LOG.info('Set -geometry flag to generate microfibrillar structure PDB file')
                return None, None
            
            return await self._handle_full_generation()
            
        except GeometryGenerationError:
            raise
        except Exception as e:
            LOG.error(f"Unexpected error in geometry generation: {str(e)}")
            traceback.print_exc()
            raise GeometryGenerationError(
                message=f"Unexpected error in geometry generation: {str(e)}",
                original_error=e,
                error_code="GEO_ERR_001", 
                context={"config": self.config.model_dump()}
            )

async def build_geometry(config: ColbuilderConfig, file_manager: Optional[FileManager] = None) -> Optional[System]:
    """
    Build geometry from configuration.
    
    Standalone function that creates a GeometryService and invokes its
    build_geometry method.
    
    Args:
        config: Configuration for geometry operations
        file_manager: Optional file manager for consistent file handling
        
    Returns:
        Generated system or None
    """
    service = GeometryService(config, file_manager)
    system, _ = await service.build_geometry()
    return system

async def build_geometry_anywhere(config: ColbuilderConfig, file_manager: Optional[FileManager] = None) -> Tuple[Path, Path]:
    """
    Main entry point for geometry generation.

    This function determines the appropriate geometry generation mode 
    based on configuration and delegates to the appropriate handler.

    Args:
        config: Configuration settings
        file_manager: Optional file manager to use for consistent file handling

    Returns:
        Tuple containing paths to output directory and final PDB file

    Raises:
        GeometryGenerationError: If geometry generation fails
    """
    try:
        if file_manager is None:
            file_manager = FileManager(config)

        if not config.pdb_file:
            raise GeometryGenerationError(
                message="PDB file not specified for geometry generation",
                error_code="GEO_ERR_005"
            )

        if not config.contact_distance and not config.crystalcontacts_file:
            raise GeometryGenerationError(
                message="Either contact_distance or crystalcontacts_file must be provided",
                error_code="GEO_ERR_001",
                context={
                    "contact_distance": config.contact_distance,
                    "crystalcontacts_file": config.crystalcontacts_file
                }
            )

        # Create a temporary directory for geometry generation
        temp_dir = file_manager.get_temp_path("geometry_gen", create_dir=True)
        LOG.info(f"{Fore.BLUE}Created temporary directory for geometry generation: {temp_dir}{Style.RESET_ALL}")

        # Change to the temporary directory
        original_dir = Path.cwd()
        os.chdir(temp_dir)

        try:
            # Copy the PDB file to the temporary directory
            pdb_path = file_manager.copy_to_directory(config.pdb_file, dest_dir=temp_dir)
            config.pdb_file = pdb_path

            # Build the geometry using the crystal builder
            crystal_builder = CrystalBuilder()
            crystal_builder.set_file_manager(file_manager)

            system = await crystal_builder.build(config)

            # Get the output path for the final PDB file
            output_pdb_path = file_manager.get_output_path(config.output, ".pdb")
            if not output_pdb_path.exists():
                raise GeometryGenerationError(
                    message=f"Expected output file not found: {output_pdb_path}",
                    error_code="GEO_ERR_001"
                )

            LOG.info(f"{Fore.BLUE}Geometry generation completed successfully. Output: {output_pdb_path}{Style.RESET_ALL}")
            return temp_dir, output_pdb_path

        finally:
            os.chdir(original_dir)

    except GeometryGenerationError as e:
        LOG.error(f"Error during geometry generation: {str(e)}")
        raise
    except Exception as e:
        LOG.error(f"Unexpected error during geometry generation: {str(e)}", exc_info=True)
        raise GeometryGenerationError(
            message=f"Failed to complete geometry generation: {str(e)}",
            original_error=e,
            error_code="GEO_ERR_001"
        )

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
    file_manager = FileManager(config)
    original_dir = Path.cwd()
    temp_dir = None
    
    try:
        temp_dir = file_manager.get_temp_path("replace_gen", create_dir=True)
        LOG.info(f"Created temporary directory for replacement: {temp_dir}")
        
        os.chdir(temp_dir)
        
        replacer = CrosslinkReplacer()
        temp_config = config.copy()
        temp_config.working_directory = temp_dir
        
        is_direct_replacement = (
            system is None and 
            config.replace_file and 
            os.path.exists(config.replace_file) and
            _is_pdb_file(config.replace_file)
        )
        
        if is_direct_replacement:
            replace_path = Path(config.replace_file)
            local_replace = temp_dir / replace_path.name
            shutil.copy2(replace_path, local_replace)
            temp_config.replace_file = local_replace
            
            LOG.info("Using direct replacement approach")
            await replacer.replace_direct(temp_config)
            
            output_prefix = temp_config.species
            if temp_config.output:
                output_prefix = temp_config.output
                
            output_pdb_path = temp_dir / f"{output_prefix}.pdb"
            LOG.info(f"PDB file written to {output_pdb_path}")
            
            final_pdb = file_manager.copy_to_output(output_pdb_path)
            
            if final_pdb.exists():
                try:
                    crystal = Crystal(pdb=str(final_pdb))
                    return System(crystal=crystal)
                except Exception as e:
                    LOG.warning(f"Could not create minimal system after direct replacement: {e}")
            
            return None
        else:
            if not system:
                raise GeometryGenerationError(
                    message="No system provided for replacement and input file is not a valid PDB",
                    error_code="GEO_ERR_004"
                )
                
            system = await replacer.replace_in_system(system, temp_config)
            
            output_prefix = temp_config.species
            if temp_config.output:
                output_prefix = temp_config.output
                
            output_pdb_path = temp_dir / f"{output_prefix}.pdb"
            system.write_pdb(pdb_out=output_prefix, fibril_length=temp_config.fibril_length, cleanup=False)
            LOG.info(f"PDB file written to {output_pdb_path}")
            
            file_manager.copy_to_output(output_pdb_path)
            
            return system
    finally:
        os.chdir(original_dir)
        
        # if not config.debug and temp_dir:
        #     file_manager.cleanup()

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