"""
This module serves as the main entry point for the sequence generation process in the ColBuilder pipeline.

It provides a high-level function, `build_sequence`, which orchestrates the generation of collagen
structures from input configurations. This includes sequence alignment, structure modeling, and
crosslink application, leveraging the `SequenceGenerator` class.

Key Features:
--------------
1. **Sequence Generation**:
   - Align input sequences with templates.
   - Generate multiple sequence alignments (MSA) and 3D protein structures.

2. **Error Handling**:
   - Raise detailed exceptions for sequence generation or system operation failures.

3. **Integration with ColBuilder**:
   - Designed to work seamlessly within the ColBuilder pipeline for collagen modeling.

Usage:
------
This module is intended to be used as part of the ColBuilder pipeline. The main function,
`build_sequence`, takes a configuration object as input and outputs the paths to the generated
MSA and PDB files.

Example:
--------
```python
from colbuilder.core.sequence.main_sequence import build_sequence
from colbuilder.core.utils.config import ColbuilderConfig

# Define configuration
config = ColbuilderConfig(
    working_directory="/path/to/working_dir",
    fasta_file="/path/to/input.fasta",
    crosslink=True,
    debug=False
)

# Generate sequence
msa_output, pdb_output = await build_sequence(config)

print(f"MSA file saved to: {msa_output}")
print(f"PDB file saved to: {pdb_output}")
```
"""

# Copyright (c) 2024, ColBuilder Development Team
# Distributed under the terms of the Apache License 2.0

from pathlib import Path
from typing import Tuple
from colorama import Fore, Style
from colbuilder.core.utils.config import ColbuilderConfig
from .sequence_generator import SequenceGenerator

from ..utils.logger import setup_logger

LOG = setup_logger(__name__)


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
    msa_output, final_pdb = await generator.generate()

    if final_pdb and final_pdb.exists():
        LOG.info(
            f"{Fore.GREEN}Sequence generation completed, output PDB: {final_pdb.absolute()}{Style.RESET_ALL}"
        )

    return msa_output, final_pdb


__all__ = ["build_sequence"]
