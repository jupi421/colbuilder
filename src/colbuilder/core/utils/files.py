"""
File management utilities for the ColBuilder system.

This module provides utilities for handling file operations, resource management, and progress 
tracking throughout the ColBuilder pipeline. It includes context managers for resource tracking, 
output suppression, PDB file manipulation, and a centralized `FileManager` class for consistent 
file handling.

Key Features:
--------------
1. **File Management**:
   - Centralized handling of temporary and permanent files.
   - Create, track, and clean up temporary files and directories.
   - Copy files to output directories with consistent logging and error handling.

2. **Context Managers**:
   - `managed_resources`: Tracks operation performance and ensures proper resource cleanup.
   - `suppress_output`: Suppresses stdout and stderr during execution.
   - Temporary file and directory context managers for safe resource handling.

3. **PDB File Manipulation**:
   - Update PDB headers to standardize structure and remove duplicates.
   - Extract and rewrite CRYST1 and ATOM records with optional custom headers.

4. **Progress Tracking**:
   - `ProgressTracker` class to monitor and log progress through multi-step operations.

5. **File Search**:
   - Locate files across multiple standard or user-defined directories.

Usage:
------
This module is designed to be used as part of the ColBuilder pipeline for managing files and 
resources. It provides utilities for temporary file handling, directory creation, and PDB 
manipulation, ensuring consistent and reliable file operations.

Example:
--------
```python
from colbuilder.core.utils.files import FileManager, managed_resources, suppress_output

# Initialize FileManager
config = ColbuilderConfig(working_directory=Path("/path/to/working_dir"))
file_manager = FileManager(config)

# Create a temporary file
with file_manager.temp_file_context("temp_file", ".txt") as temp_file:
    with open(temp_file, "w") as f:
        f.write("Temporary content")

# Suppress output during noisy operations
with suppress_output():
    print("This will not appear in the console")

# Track operation performance
with managed_resources("Example Operation"):
    # Perform some operation
    pass
```
"""

import contextlib
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from pathlib import Path
from typing import Generator, Optional, Protocol, TypeVar, Any, Set, Dict, List
import io
import time
import shutil
import os
from dataclasses import dataclass

from colbuilder.core.utils.exceptions import (
    SequenceGenerationError, 
    SystemError
)
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

OperationType = TypeVar('OperationType') 

@dataclass
class OperationContext:
    """
    Base context for all operations.
    
    Attributes:
        config: Configuration for the operation
        working_dir: Working directory path
    """
    config: ColbuilderConfig
    working_dir: Path


class Operation(Protocol[OperationType]):
    """
    Base protocol for all operations.
    
    This protocol defines the common interface that all operations must implement,
    allowing for consistent execution patterns throughout the system.
    """
    
    async def execute(self, input_data: Optional[Any] = None) -> OperationType:
        """
        Execute the operation.
        
        Args:
            input_data: Optional input data for the operation
            
        Returns:
            The operation result of type OperationType
        """
        pass


@contextmanager
def managed_resources(operation_name: str) -> Generator[None, None, None]:
    """
    Context manager for tracking operation performance and managing resources.
    
    This context manager handles logging of operation start/end times and ensures
    proper resource cleanup even if exceptions occur during operation execution.
    
    Args:
        operation_name: Name of the operation being performed
        
    Yields:
        None
        
    Raises:
        SystemError: If resource management fails
    """
    start_time = time.perf_counter()
    try:
        LOG.debug(f"Starting operation: {operation_name}")
        yield
    except Exception as e:
        LOG.error(f"Error in operation {operation_name}: {str(e)}")
        raise
    finally:
        try:
            duration = time.perf_counter() - start_time
            LOG.debug(f"{operation_name} completed in {duration:.2f} seconds")
        except Exception as e:
            raise SystemError(
                message="Failed to finalize resource management",
                original_error=e,
                error_code="SYS_ERR_001",
                context={
                    "operation": operation_name,
                    "duration": time.perf_counter() - start_time
                }
            )


@contextmanager
def suppress_output() -> Generator[None, None, None]:
    """
    Context manager to suppress stdout and stderr output.
    
    This utility captures and discards all standard output and error streams
    during its execution scope, which is useful when calling noisy external
    libraries or tools where their console output is not relevant.
    
    Yields:
        None
        
    Raises:
        SystemError: If there's an error managing output streams
    """
    try:
        with io.StringIO() as stdout_buf, io.StringIO() as stderr_buf:
            with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
                yield
    except Exception as e:
        raise SystemError(
            message="Failed to manage output streams",
            original_error=e,
            error_code="SYS_ERR_001",
            context={"action": "suppress_output"}
        )


class ProgressTracker:
    """
    Tracks progress of multi-step operations.
    
    This class provides a simple way to track and report progress through
    a sequence of steps, maintaining state about how far along the process is.
    """
    
    def __init__(self, total_steps: int) -> None:
        """
        Initialize the progress tracker.
        
        Args:
            total_steps: Total number of steps in the operation
        """
        self.total_steps: int = total_steps
        self.current_step: int = 0
    
    def update(self, message: str) -> None:
        """
        Update progress with a step completion message.
        
        Increments the step counter and logs the progress message with
        the current step number and total steps.
        
        Args:
            message: Progress message to log
        """
        self.current_step += 1
        LOG.info(f"Step {self.current_step}/{self.total_steps}: {message}")


def update_pdb_header(pdb_file: Path, first_line: str = '') -> None:
    """
    Update PDB file header, removing duplicates and unnecessary lines.
    
    This function reads a PDB file, extracts CRYST1 and ATOM records,
    and writes them back with an optional custom first line. It helps
    standardize PDB files and ensure they have the correct header structure.
    
    Args:
        pdb_file: Path to PDB file
        first_line: Optional first line to add to the file
    
    Raises:
        SequenceGenerationError: If update fails
    """
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        cryst_lines = []
        atom_lines = []
        in_atom_section = False
        
        for line in lines:
            if line.startswith('CRYST1'):
                cryst_lines.append(line)
            elif line.startswith('ATOM'):
                in_atom_section = True
                atom_lines.append(line)
            elif in_atom_section:
                atom_lines.append(line)
        
        output_lines = []
        
        if first_line and not first_line.isspace() and len(first_line) > 0:
            output_lines.append(first_line + '\n')
        
        if cryst_lines:
            output_lines.append(cryst_lines[0])
            
        output_lines.extend(atom_lines)
        
        with open(pdb_file, 'w') as f:
            f.writelines(output_lines)
            
    except Exception as e:
        LOG.error(f"Error updating PDB header: {str(e)}")
        raise SequenceGenerationError(
            "Failed to update PDB header",
            error_code="SEQ_ERR_004",
            context={
                "pdb_file": str(pdb_file),
                "error": str(e)
            }
        )


class FileManager:
    """
    Centralized utility for file management in Colbuilder.
    
    This class provides consistent methods for managing temporary and permanent files,
    ensuring proper cleanup and respecting debug settings. It handles file paths,
    temporary directories, copying files, and directory creation with consistent
    logging and error handling.
    
    Attributes:
        config: Colbuilder configuration
        temp_files: Set of temporary files to track
        temp_dirs: Set of temporary directories to track
    """
    
    def __init__(self, config: ColbuilderConfig) -> None:
        """
        Initialize with configuration settings.
        
        Args:
            config: Colbuilder configuration
        """
        self.config: ColbuilderConfig = config
        self.project_root = config.PROJECT_ROOT  # Use the project root from config
        self.working_dir = Path(config.working_directory).resolve()
        self.temp_files: Set[Path] = set()
        self.temp_dirs: Set[Path] = set()
        
    def get_temp_path(self, basename: str, suffix: Optional[str] = None, create_dir: bool = False) -> Path:
        """
        Get a path for a temporary file or directory.
        
        Creates a path for a temporary file or directory, optionally creating
        the directory structure. The path is tracked for later cleanup.
        
        Args:
            basename: Base name for the file or directory
            suffix: Optional suffix to add (e.g., file extension)
            create_dir: Whether to create a directory at this path
            
        Returns:
            Path object for the temporary file or directory
        """
        if self.config.debug:
            base_dir = self.config.working_directory / "temp_working_dir"
        else:
            base_dir = self.config.working_directory / ".tmp"
            
        base_dir.mkdir(exist_ok=True, parents=True)
        self.temp_dirs.add(base_dir)
        
        full_name = f"{basename}{suffix if suffix else ''}"
        full_path = base_dir / full_name
        
        if create_dir:
            full_path.mkdir(exist_ok=True, parents=True)
            self.temp_dirs.add(full_path)
        else:
            self.temp_files.add(full_path)
            
        return full_path
    
    def get_output_path(self, basename: str, suffix: Optional[str] = None, mkdir: bool = False) -> Path:
        """
        Get a path for a permanent output file or directory.
        
        Creates a path for an output file or directory in the working directory,
        optionally creating the directory structure.
        
        Args:
            basename: Base name for the file or directory
            suffix: Optional suffix to add (e.g., file extension)
            mkdir: Whether to create a directory at this path
            
        Returns:
            Path object for the output file or directory
        """
        full_name = f"{basename}{suffix if suffix else ''}"
        full_path = self.config.working_directory / full_name
        
        if mkdir:
            full_path.mkdir(exist_ok=True, parents=True)
            
        return full_path
    
    def cleanup(self, force: bool = False) -> None:
        """
        Clean up temporary files and directories.
        
        Removes all tracked temporary files and directories, unless
        in debug mode and not forced.
        
        Args:
            force: Whether to clean up even in debug mode
        """
        if self.config.debug and not force:
            LOG.debug("Skipping cleanup due to debug mode")
            return
            
        for file_path in self.temp_files:
            try:
                if file_path.exists():
                    file_path.unlink()
                    LOG.debug(f"Removed temporary file: {file_path}")
            except Exception as e:
                LOG.warning(f"Failed to remove temporary file {file_path}: {str(e)}")
                
        for dir_path in sorted(self.temp_dirs, key=lambda p: -len(str(p))):
            try:
                if dir_path.exists():
                    shutil.rmtree(dir_path)
                    LOG.debug(f"Removed temporary directory: {dir_path}")
            except Exception as e:
                LOG.warning(f"Failed to remove temporary directory {dir_path}: {str(e)}")
                
        self.temp_files.clear()
        self.temp_dirs.clear()
        
    def copy_to_output(self, source: Path, dest_name: Optional[str] = None) -> Path:
        """
        Copy a file to the output directory.
        
        Copies a source file to the working directory with the given
        destination name, or the original name if not specified.
        
        Args:
            source: Source file path
            dest_name: Optional destination name (uses source name if not provided)
            
        Returns:
            Path to the destination file
            
        Raises:
            Various exceptions from shutil.copy2 if copy fails
        """
        dest_name = dest_name or source.name
        dest_path = self.config.working_directory / dest_name
        
        try:
            shutil.copy2(source, dest_path)
            LOG.debug(f"Copied {source} to {dest_path}")
            return dest_path
        except Exception as e:
            LOG.error(f"Failed to copy {source} to {dest_path}: {str(e)}")
            raise
            
    def ensure_dir(self, dirname: str) -> Path:
        """
        Ensure a directory exists in the working directory.
        
        Creates the directory if it doesn't already exist.
        
        Args:
            dirname: Name of the directory
            
        Returns:
            Path to the directory
        """
        dir_path = self.config.working_directory / dirname
        dir_path.mkdir(exist_ok=True, parents=True)
        return dir_path

    @contextmanager
    def temp_file_context(self, basename: str, suffix: Optional[str] = None) -> Generator[Path, None, None]:
        """
        Context manager for temporary file operations.
        
        Creates a temporary file path, yields it for use, and ensures
        it gets tracked for cleanup later.
        
        Args:
            basename: Base name for the file
            suffix: Optional suffix/extension
            
        Yields:
            Path to the temporary file
        """
        temp_path = self.get_temp_path(basename, suffix)
        try:
            yield temp_path
        finally:
            self.temp_files.add(temp_path)
            
    @contextmanager
    def temp_dir_context(self, dirname: str) -> Generator[Path, None, None]:
        """
        Context manager for temporary directory operations.
        
        Creates a temporary directory, yields it for use, and ensures
        it gets tracked for cleanup later.
        
        Args:
            dirname: Name for the directory
            
        Yields:
            Path to the temporary directory
        """
        temp_dir = self.get_temp_path(dirname, create_dir=True)
        try:
            yield temp_dir
        finally:
            self.temp_dirs.add(temp_dir)
            
    def find_file(self, filename_or_path: str, search_paths: Optional[List[Path]] = None) -> Optional[Path]:
        """
        Find a file by searching in multiple locations.
        
        Searches for a file in the provided search paths, standard Colbuilder locations,
        and relative to the project root.
        
        Args:
            filename_or_path: Name of the file or relative path to find
            search_paths: Optional list of paths to search (added to standard locations)
            
        Returns:
            Path to the file if found, None otherwise
        """
        # Convert string to Path if necessary
        filepath = Path(filename_or_path)
        
        # First, check if it's an absolute path
        if filepath.is_absolute() and filepath.exists():
            return filepath
            
        # Initialize standard search paths if none provided
        if search_paths is None:
            search_paths = []
            
        # Add standard search locations
        search_paths.extend([
            self.working_dir,                 # User's current working directory
            Path.cwd(),                       # Python process working directory
            self.project_root,                # Project root from config
            self.project_root / "data",       # Project data directory
            self.project_root / "colbuilder" / "data",  # Package data directory
            self.config.DATA_DIR,             # Data directory from config
            self.config.HOMOLOGY_LIB_DIR,     # Homology lib directory from config
            self.config.FORCE_FIELD_DIR       # Force field directory from config
        ])
        
        # Filter out None values
        search_paths = [p for p in search_paths if p is not None]
        
        # Try finding the file in each search path
        for base_path in search_paths:
            # Try with the filepath as provided
            full_path = base_path / filepath
            if full_path.exists():
                return full_path
                
            # If filepath has multiple parts, also try with just the filename
            if len(filepath.parts) > 1:
                filename_only = filepath.name
                name_only_path = base_path / filename_only
                if name_only_path.exists():
                    return name_only_path
                
        # File not found
        return None

    def copy_to_directory(self, source: Path, dest_dir: Optional[Path] = None, dest_name: Optional[str] = None) -> Path:
        """
        Copy a file to a specified directory.

        Args:
            source (Path): Source file path.
            dest_dir (Optional[Path]): Destination directory (defaults to the working directory).
            dest_name (Optional[str]): Optional destination name (uses source name if not provided).

        Returns:
            Path: Path to the copied file in the destination directory.

        Raises:
            FileNotFoundError: If the source file does not exist.
            Exception: If the copy operation fails.
        """
        if not source.exists():
            raise FileNotFoundError(f"Source file not found: {source}")

        dest_dir = dest_dir or self.config.working_directory
        dest_dir.mkdir(parents=True, exist_ok=True)

        dest_name = dest_name or source.name
        dest_path = dest_dir / dest_name

        try:
            shutil.copy2(source, dest_path)
            LOG.debug(f"Copied {source} to {dest_path}")
            return dest_path
        except Exception as e:
            LOG.error(f"Failed to copy {source} to {dest_path}: {str(e)}")
            raise

    def get_type_dir(self, model_type: str) -> Path:
        """
        Get or create the directory for a specific model type.

        Args:
            model_type (str): The type of the model (e.g., "D", "NC").

        Returns:
            Path: Path to the directory for the specified model type.

        Raises:
            FileNotFoundError: If the directory cannot be created or accessed.
        """
        try:
            # Define the base directory for model types
            type_dir = self.config.working_directory / "types" / model_type

            # Create the directory if it doesn't exist
            type_dir.mkdir(parents=True, exist_ok=True)

            return type_dir
        except Exception as e:
            raise FileNotFoundError(f"Failed to locate or create type directory for model type: {model_type}") from e

    def get_temp_dir(self, dirname: str) -> Path:
        """
        Get a temporary directory path and create the directory.
        
        This is a convenience wrapper around get_temp_path with create_dir=True.
        
        Args:
            dirname: Name for the directory
            
        Returns:
            Path to the created temporary directory
        """
        temp_dir = self.get_temp_path(dirname, create_dir=True)
        LOG.info(f"Created temporary directory: {temp_dir}")
        return temp_dir