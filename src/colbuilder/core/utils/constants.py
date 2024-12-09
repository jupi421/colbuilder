# constants.py
from typing import Final
from pathlib import Path

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