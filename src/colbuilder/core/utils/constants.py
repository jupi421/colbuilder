"""
Constants for the ColBuilder Pipeline

This module defines constants used throughout the ColBuilder pipeline for optimization, file 
operations, and PDB formatting. These constants provide a centralized and consistent reference 
to ensure uniform behavior across the system.

Key Features:
--------------
1. **Optimization Constants**:
   - `MAX_OPTIMIZATION_ATTEMPTS`: Maximum number of optimization attempts for crosslinks.
   - `MAX_TRIVALENT_DISTANCE`: Maximum allowable distance for trivalent crosslinks (in Ångstroms).
   - `MAX_DIVALENT_DISTANCE`: Maximum allowable distance for divalent crosslinks (in Ångstroms).
   - `CRITICAL_DISTANCE_THRESHOLD`: Distance threshold beyond which optimization is considered critical.

2. **File Operation Constants**:
   - `TEMP_FILE_SUFFIX`: Suffix for temporary files.
   - `DISORIENTED_SUFFIX`: Suffix for disoriented files.
   - `PDB_EXTENSION`: File extension for PDB files.
   - `FASTA_EXTENSION`: File extension for FASTA files.

3. **PDB Format Constants**:
   - `DEFAULT_PDB_HEADER`: Default header for PDB files.
   - `ATOM_RECORD_LENGTH`: Standard length of ATOM records in PDB files.

Usage:
------
This module is intended to be imported wherever these constants are required to ensure consistency 
and avoid hardcoding values.

Example:
--------
```python
from colbuilder.core.utils.constants import MAX_OPTIMIZATION_ATTEMPTS, PDB_EXTENSION

# Use constants in file operations
temp_file = f"structure{TEMP_FILE_SUFFIX}{PDB_EXTENSION}"
print(f"Temporary file: {temp_file}")

# Use constants in optimization logic
if attempts > MAX_OPTIMIZATION_ATTEMPTS:
    raise ValueError("Exceeded maximum optimization attempts.")
```
"""
from typing import Final

# Optimization constants
MAX_OPTIMIZATION_ATTEMPTS: Final = 3
MAX_TRIVALENT_DISTANCE: Final = 7.0
MAX_DIVALENT_DISTANCE: Final = 5.0
CRITICAL_DISTANCE_THRESHOLD: Final = 9.0

# File operation constants
TEMP_FILE_SUFFIX: Final = "_temp"
DISORIENTED_SUFFIX: Final = "_disoriented"
PDB_EXTENSION: Final = ".pdb"
FASTA_EXTENSION: Final = ".fasta"

# PDB format constants
DEFAULT_PDB_HEADER: Final = "CRYST1   39.970   26.950  677.900  89.24  94.59 105.58 P 1           2"
ATOM_RECORD_LENGTH: Final = 80