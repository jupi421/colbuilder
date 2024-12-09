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

def update_pdb_header(file_path: Path, header: str) -> None:
    """
    Add a new header to a PDB file and remove REMARK lines.
    
    Args:
        file_path: Path to the PDB file
        header: New header to add
        
    Raises:
        SequenceGenerationError: If file operations fail
    """
    try:
        with managed_resources("PDB Header Update"):  
            with open(file_path, "r") as f:
                lines = f.readlines()
            
            lines.insert(0, header + '\n')
            lines = [line for line in lines if not line.startswith("REMARK")]
            
            with open(file_path, "w") as f:
                f.writelines(lines)
                
    except IOError as e:
        raise SequenceGenerationError(
            message="Failed to update PDB header",
            original_error=e,
            error_code="SEQ_ERR_001",
            context={
                "file_path": str(file_path),
                "operation": "update_header"
            }
        )