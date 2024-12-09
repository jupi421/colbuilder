from pathlib import Path
from typing import Tuple
from colbuilder.core.utils.config import ColbuilderConfig
from .sequence_generator import SequenceGenerator

async def build_sequence(config: ColbuilderConfig) -> Tuple[Path, Path]:
    """
    Entry point for sequence generation process.
    
    Args:
        config: Configuration object
        
    Returns:
        Tuple of (msa_output_path, final_pdb_path)
        
    Raises:
        SequenceGenerationError: If sequence generation fails
        SystemError: If system operations fail
    """
    generator = SequenceGenerator(config)
    return await generator.generate()

__all__ = ['build_sequence']