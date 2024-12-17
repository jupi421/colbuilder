# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
import numpy as np
from typing import List, Optional, Union
from pathlib import Path

class Crosslink:
    """
    Container class storing information about each model's crosslink.

    Attributes:
        model_id (Optional[Union[int, float]]): ID of the model containing this crosslink.
        resid (str): Residue ID.
        resname (str): Residue name.
        chain (str): Chain identifier.
        position (np.ndarray): 3D coordinates of the crosslink.
        type (str): Type of the crosslink ('T' or 'D').
        state (str): State of the crosslink (default is 'none').
    """

    def __init__(self, resid: str, resname: str, chain: str, position: List[float], 
                 type: str, model_id: Optional[Union[int, float]] = None):
        self.model_id = model_id
        self.resid = resid
        self.resname = resname
        self.chain = chain
        self.position = np.array(position)
        self.type = type
        self.state = 'none'

    def set_transform(self, transform: List[float], model_id: Union[int, float]) -> None:
        """
        Sets model ID and translates crosslink according to transformation matrix.

        Args:
            transform (List[float]): Translation vector.
            model_id (Union[int, float]): ID of the model.
        """
        self.model_id = model_id
        self.position = np.add(transform, self.position)

    def __repr__(self) -> str:
        """
        String representation of the Crosslink object.

        Returns:
            str: String representation of the Crosslink.
        """
        return (f"Crosslink(model_id={self.model_id}, resid={self.resid}, "
                f"resname={self.resname}, chain={self.chain}, "
                f"position={self.position}, type={self.type}, state={self.state})")


def read_crosslink(pdb_file: Union[str, Path]) -> List[Crosslink]:
    """
    Reads the crosslink information from a PDB file.

    Args:
        pdb_file (Union[str, Path]): Path to the PDB file.

    Returns:
        List[Crosslink]: List of Crosslink objects.
    """
    crosslinks = []
    pdb_path = Path(pdb_file).with_suffix('.pdb')
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if ((line[17:20] == ('LYX' or 'LXY' or 'LYY' or 'LXX') and line[13:16] == ('C13' or 'C12')) or
                (line[17:20] == ('LY3' or 'LX3' or 'L3Y' or 'L2Y' or 'L3X' or 'L2X') and line[13:15] == 'CG') or
                (line[17:20] == ('LY2' or 'LX2') and line[13:15] == 'CB')):
                crosslinks.append(Crosslink(
                    resid=line[22:26].strip(),
                    resname=line[17:20],
                    chain=line[21],
                    position=[float(line[29:38]), float(line[38:46]), float(line[46:56])],
                    type='T'
                ))
            elif ((line[17:20] == ('L4Y' or 'L4X' or 'LY4' or 'LX4') and line[13:15] == 'CE') or
                  (line[17:20] == ('L5Y' or 'L5X' or 'LY5' or 'LX5') and line[13:15] == 'NZ')):
                crosslinks.append(Crosslink(
                    resid=line[22:26].strip(),
                    resname=line[17:20],
                    chain=line[21],
                    position=[float(line[29:38]), float(line[38:46]), float(line[46:56])],
                    type='D'
                ))
            elif ((line[17:20] == ('LGX' or 'LPS') and line[13:15] == 'CE') or
                  (line[17:20] == ('AGS' or 'APD') and line[13:15] == 'NZ')):
                crosslinks.append(Crosslink(
                    resid=line[22:26].strip(),
                    resname=line[17:20],
                    chain=line[21],
                    position=[float(line[29:38]), float(line[38:46]), float(line[46:56])],
                    type='D'
                ))
   
    return crosslinks
