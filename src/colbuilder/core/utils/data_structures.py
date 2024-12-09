# data_structures.py
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import numpy as np

@dataclass(frozen=True)
class CrosslinkPosition:
    """
    Represents a position in a crosslink, including residue and atom information.
    
    Attributes:
        residue_number (int): The residue sequence number
        chain_id (str): Chain identifier (A, B, C)
        residue_type (str): Three-letter residue code
        atom_name (str): Atom identifier in the residue
        
    Raises:
        ValueError: If chain_id is not A, B, or C
    """
    residue_number: int
    chain_id: str
    residue_type: str
    atom_name: str
    
    def __post_init__(self):
        if self.chain_id not in {'A', 'B', 'C'}:
            raise ValueError(f"Chain ID must be A, B, or C, got {self.chain_id}")
        
    @property
    def position_str(self) -> str:
        """Returns the position in format 'number.chain'"""
        return f"{self.residue_number}.{self.chain_id}"

@dataclass(frozen=True)
class CrosslinkPair:
    """
    Represents a crosslink between two or three positions.
    
    Attributes:
        position1 (CrosslinkPosition): First position in crosslink
        position2 (CrosslinkPosition): Second position in crosslink
        position3 (Optional[CrosslinkPosition]): Optional third position for trivalent crosslinks
        terminal_type (str): Terminal type (N or C)
    """
    position1: CrosslinkPosition
    position2: CrosslinkPosition
    position3: Optional[CrosslinkPosition] = None
    terminal_type: str = field(default="N")
    
    def __post_init__(self):
        if self.terminal_type not in {'N', 'C'}:
            raise ValueError(f"Terminal type must be N or C, got {self.terminal_type}")
            
    @property
    def is_trivalent(self) -> bool:
        """Returns True if this is trivalent crosslink."""
        return self.position3 is not None

@dataclass
class OptimizationState:
    """
    Tracks the state of crosslink optimization.
    
    Attributes:
        current_distance: Current distance between crosslinked atoms
        attempt_number: Current optimization attempt number
        best_distance: Best distance achieved so far
        coordinates: Numpy array of current atomic coordinates
    """
    current_distance: float = float('inf')
    attempt_number: int = 0
    best_distance: float = float('inf')
    coordinates: Optional[np.ndarray] = None
    optimization_history: List[Dict[str, Any]] = field(default_factory=list)
    
    def update(self, distance: float, coords: Optional[np.ndarray] = None) -> None:
        """
        Update optimization state with new values.
        
        Args:
            distance: Current distance between crosslinked atoms
            coords: Optional coordinates array. If None, coordinates won't be updated.
        """
        self.current_distance = distance
        if distance < self.best_distance:
            self.best_distance = distance
            if coords is not None:
                self.coordinates = coords.copy()
        self.optimization_history.append({
            'attempt': self.attempt_number,
            'distance': distance
        })
        
    def increment_attempt(self) -> None:
        """Increment the attempt counter."""
        self.attempt_number += 1