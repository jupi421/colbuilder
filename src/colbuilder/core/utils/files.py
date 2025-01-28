# files.py
import contextlib
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from pathlib import Path
from typing import Generator, Optional, Protocol, TypeVar, Any
import io
import time
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
    """Base context for all operations."""
    config: ColbuilderConfig
    working_dir: Path

class Operation(Protocol[OperationType]):
    """Base protocol for all operations."""
    async def execute(self, input_data: Optional[Any] = None) -> OperationType:
        """Execute the operation."""
        pass

@contextmanager
def managed_resources(operation_name: str) -> Generator[None, None, None]:
    """
    Context manager for tracking operation performance and managing resources.
    
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
    Context manager to suppress stdout and stderr.
    
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
    """Tracks progress of multi-step operations."""
    
    def __init__(self, total_steps: int) -> None:
        """
        Initialize progress tracker.
        
        Args:
            total_steps: Total number of steps in the operation
        """
        self.total_steps = total_steps
        self.current_step = 0
    
    def update(self, message: str) -> None:
        """
        Update progress with a message.
        
        Args:
            message: Progress message to log
        """
        self.current_step += 1
        LOG.info(f"Step {self.current_step}/{self.total_steps}: {message}")

def update_pdb_header(pdb_file: Path, first_line: str = '') -> None:
    """
    Update PDB file header, removing duplicates and unnecessary lines.
    
    Args:
        pdb_file: Path to PDB file
        first_line: Optional first line to add to the file
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