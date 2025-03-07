# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import numpy as np
from sklearn.metrics import pairwise_distances as pdist
from typing import List, Dict, Any, Optional, Tuple, Union
import os

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class Crosslink:
    """
    Setup crosslink topology for collagen models.
    
    Identifies crosslink sites in merged PDB files and generates
    the necessary bonded parameters (bonds, angles, dihedrals) for crosslinks
    in Martini coarse-grained models.
    """
    
    def __init__(self, cnt_model: Optional[int] = None):
        """
        Initialize the Crosslink object.
        
        Parameters
        ----------
        cnt_model : Optional[int]
            The model counter for file naming
        """
        self.file = f"{int(cnt_model)}.merge.pdb" if cnt_model is not None else None
        self.crosslink_coords: List[List[float]] = []
        self.crosslink_pdb: List[List[str]] = []
        self.crosslink_neighbors: List[Any] = []
        self.crosslink_connect: List[List[List[str]]] = []
        self.crosslink_bonded: Dict[str, List[List[Any]]] = {'bonds': [], 'angles': [], 'dihedrals': []}
        
        # Parameters for the divalent HLKNL-crosslink
        self.dly45 = '0.415'  # Equilibrium distance for L4Y-L5Y bond
        self.kly45 = '7000'   # Force constant for L4Y-L5Y bond
        self.al45y_1 = '140'  # L4Y-L5Y link SC1-SC2-SC1(L4Y) - equilibrium angle
        self.al45y_2 = '140'  # L4Y-L5Y link SC2(L5Y)-SC1-BB - equilibrium angle
        self.k_angle = '153'  # Force constant for all angles
        
        # Parameters for the trivalent PYD-crosslink
        self.klyxly2 = '9000'   # Force constant for LYX-LY2 bond
        self.klyxly3 = '12000'  # Force constant for LYX-LY3 bond
        self.dlyxly2 = '0.290'  # Equilibrium distance for LYX-LY2 bond
        self.dlyxly3 = '0.230'  # Equilibrium distance for LYX-LY3 bond
        
        self.al2yx_1 = '100'  # LY2-LYX link TP1q-TC6q-TC4 - equilibrium angle
        self.al2yx_2 = '60'   # LY2-LYX link TQ2p-TP1q-TC4 - equilibrium angle
        self.al2yx_3 = '130'  # LY2-LYX link TP1q-TC4-SP2 - equilibrium angle
        self.al3yx_1 = '100'  # LY3-LYX link TP1q-TC6q-TC4 - equilibrium angle
        self.al3yx_2 = '110'  # LY3-LYX link TQ2p-TC6q-TC4 - equilibrium angle
        self.al3yx_3 = '130'  # LY3-LYX link TC6q-TC4-SP2 - equilibrium angle
        
        LOG.debug(f"Initialized Crosslink for model counter {cnt_model}")

    def get_crosslink_coords(self, cnt_model: Optional[int] = None) -> List[List[float]]:
        """
        Identify crosslinks from PDB and extract their coordinates.
        
        Reads the PDB file and extracts coordinates of atoms that
        are part of crosslinks (LYX, LY2, LY3, L4Y, L5Y residues and specific atoms).
        
        Parameters
        ----------
        cnt_model : Optional[int]
            The model counter for file naming if different from initialization
            
        Returns
        -------
        List[List[float]]
            List of coordinates for crosslink sites
        """
        if cnt_model is None:
            file = self.file
        else:
            file = f"{int(cnt_model)}.merge.pdb"
            
        LOG.debug(f"Reading crosslink coordinates from {file}")
        
        if not os.path.exists(file):
            LOG.error(f"PDB file not found: {file}")
            return self.crosslink_coords
            
        self.crosslink_coords = []
        self.crosslink_pdb = []
        it_pdb = 0
        
        try:
            with open(file, 'r') as f:
                for line in f:
                    if line[0:4] == 'ATOM':
                        it_pdb += 1
                        
                        if ((line[17:20] == 'LYX' and line[12:15] == 'SC4') or
                            (line[17:20] == 'LY2' and line[12:15] == 'SC1') or
                            (line[17:20] == 'LY3' and line[12:15] == 'SC1') or
                            (line[17:20] == 'LYX' and line[12:15] == 'SC5')):
                            
                            coords = [
                                float(line[29:38]),
                                float(line[38:46]),
                                float(line[46:56])
                            ]
                            
                            self.crosslink_pdb.append([
                                str(it_pdb),
                                line[17:20],  # residue
                                line[12:15],  # atom
                                line[21:26],  # chain/residue number
                                coords[0], coords[1], coords[2]
                            ])
                            self.crosslink_coords.append(coords)
                            
                        elif ((line[17:20] == 'L4Y' and line[12:15] == 'SC1') or
                              (line[17:20] == 'L5Y' and line[12:15] == 'SC2')):
                            
                            coords = [
                                float(line[29:38]),
                                float(line[38:46]),
                                float(line[46:56])
                            ]
                            
                            self.crosslink_pdb.append([
                                str(it_pdb),
                                line[17:20],  # residue
                                line[12:15],  # atom
                                line[21:26],  # chain/residue number
                                coords[0], coords[1], coords[2]
                            ])
                            self.crosslink_coords.append(coords)
            
            LOG.debug(f"Found {len(self.crosslink_coords)} crosslink sites in {file}")
                
        except Exception as e:
            LOG.error(f"Error reading PDB file {file}: {str(e)}")
            
        return self.crosslink_coords

    def get_crosslink_connect(self, cnt_model: Optional[int] = None) -> List[List[List[str]]]:
        """
        Get nearest crosslinks to determine connections between them.
        
        Calculates pairwise distances between crosslink sites and
        identifies the closest sites to each other.
        
        Parameters
        ----------
        cnt_model : Optional[int]
            The model counter for file naming if different from initialization
            
        Returns
        -------
        List[List[List[str]]]
            List of connected crosslink sites
        """
        LOG.debug(f"Finding crosslink connections for model {cnt_model}")
        
        self.get_crosslink_coords(cnt_model=cnt_model)
        
        if not self.crosslink_coords:
            LOG.warning("No crosslink coordinates found")
            return []
        
        try:
            pairs = pdist(self.crosslink_coords)
            out = []  # Indices that have been processed
            self.crosslink_connect = []
            
            for p in pairs:
                tmp = []
                for k in np.argsort(p)[0:4]:
                    if k not in out:
                        tmp.append(self.crosslink_pdb[k])
                        out.append(k)
                        
                if tmp:
                    self.crosslink_connect.append(tmp)
                    
            LOG.debug(f"Found {len(self.crosslink_connect)} potential crosslink connections")
                
        except Exception as e:
            LOG.error(f"Error finding crosslink connections: {str(e)}")
            
        return self.crosslink_connect
    
    def set_crosslink_bonded(self, cnt_model: Optional[int] = None, 
                           crosslink_connect: Optional[List[List[List[str]]]] = None) -> Dict[str, List[List[Any]]]:
        """
        Setup topology for crosslink bonded parameters: bonds, angles, and dihedrals.
        
        Identifies different types of crosslinks and sets up the
        appropriate bonded parameters for each.
        
        Parameters
        ----------
        cnt_model : Optional[int]
            The model counter for file naming if different from initialization
        crosslink_connect : Optional[List[List[List[str]]]]
            List of connected crosslink sites (optional)
            
        Returns
        -------
        Dict[str, List[List[Any]]]
            Dictionary of bonded parameters for crosslinks
        """
        LOG.debug(f"Setting up crosslink bonded parameters for model {cnt_model}")
        
        self.crosslink_bonded = {'bonds': [], 'angles': [], 'dihedrals': []}
        
        if crosslink_connect is None:
            crosslink_connect = self.get_crosslink_connect(cnt_model=cnt_model)
            
        if not crosslink_connect:
            LOG.debug("No crosslink connections found, returning empty parameters")
            return self.crosslink_bonded
            
        try:
            connections_found = 0
            for c in crosslink_connect:
                for i, clx in enumerate(c):
                    for j, cly in enumerate(c):
                        if i == j:
                            continue
                            
                        try:
                            dist = np.linalg.norm(np.array(clx[-3:]) - np.array(cly[-3:]))
                        except (ValueError, TypeError) as e:
                            LOG.warning(f"Error calculating distance: {str(e)}")
                            continue
                            
                        if dist < 10.0:
                            if clx[1] == 'LYX' and clx[2] == 'SC4' and cly[1] == 'LY2':
                                self._add_lyx_ly2_bonds(clx, cly)
                                connections_found += 1
                                LOG.debug(f"Added LYX-LY2 bonds between {clx[0]} and {cly[0]}")
                                
                            elif clx[1] == 'LYX' and clx[2] == 'SC5' and cly[1] == 'LY3':
                                self._add_lyx_ly3_bonds(clx, cly)
                                connections_found += 1
                                LOG.debug(f"Added LYX-LY3 bonds between {clx[0]} and {cly[0]}")
                                
                            if clx[1] == 'L4Y' and clx[2] == 'SC1' and cly[1] == 'L5Y' and cly[2] == 'SC2':
                                self._add_l4y_l5y_bonds(clx, cly)
                                connections_found += 1
                                LOG.debug(f"Added L4Y-L5Y bonds between {clx[0]} and {cly[0]}")
                    
            LOG.debug(f"Created {len(self.crosslink_bonded['bonds'])} bonds and "
                    f"{len(self.crosslink_bonded['angles'])} angles from "
                    f"{connections_found} crosslink connections")
                
        except Exception as e:
            LOG.error(f"Error setting up crosslink bonded parameters: {str(e)}")
            
        return self.crosslink_bonded
    
    def _add_lyx_ly2_bonds(self, clx: List[Any], cly: List[Any]) -> None:
        """
        Add bonded parameters for LYX-LY2 crosslink.
        
        Parameters
        ----------
        clx : List[Any]
            First crosslink atom data
        cly : List[Any]
            Second crosslink atom data
        """
        self.crosslink_bonded['bonds'].append([
            clx[0], cly[0], '1', self.dlyxly2, f"{self.klyxly2}\n"
        ])
        
        self.crosslink_bonded['angles'].append([
            str(int(clx[0])+1), clx[0], cly[0], '1', self.al2yx_1, f"{self.k_angle}\n"
        ])
        self.crosslink_bonded['angles'].append([
            str(int(clx[0])-1), clx[0], cly[0], '1', self.al2yx_2, f"{self.k_angle}\n"
        ])
        self.crosslink_bonded['angles'].append([
            clx[0], cly[0], str(int(cly[0])-1), '1', self.al2yx_3, f"{self.k_angle}\n"
        ])
    
    def _add_lyx_ly3_bonds(self, clx: List[Any], cly: List[Any]) -> None:
        """
        Add bonded parameters for LYX-LY3 crosslink.
        
        Parameters
        ----------
        clx : List[Any]
            First crosslink atom data
        cly : List[Any]
            Second crosslink atom data
        """
        self.crosslink_bonded['bonds'].append([
            clx[0], cly[0], '1', self.dlyxly3, f"{self.klyxly3}\n"
        ])
        
        self.crosslink_bonded['angles'].append([
            str(int(clx[0])-1), clx[0], cly[0], '1', self.al3yx_1, f"{self.k_angle}\n"
        ])
        self.crosslink_bonded['angles'].append([
            str(int(clx[0])-2), clx[0], cly[0], '1', self.al3yx_2, f"{self.k_angle}\n"
        ])
        self.crosslink_bonded['angles'].append([
            clx[0], cly[0], str(int(cly[0])-1), '1', self.al3yx_3, f"{self.k_angle}\n"
        ])
    
    def _add_l4y_l5y_bonds(self, clx: List[Any], cly: List[Any]) -> None:
        """
        Add bonded parameters for L4Y-L5Y crosslink.
        
        Parameters
        ----------
        clx : List[Any]
            First crosslink atom data
        cly : List[Any]
            Second crosslink atom data
        """
        self.crosslink_bonded['bonds'].append([
            clx[0], cly[0], '1', self.dly45, f"{self.kly45}\n"
        ])
        
        self.crosslink_bonded['angles'].append([
            clx[0], cly[0], str(int(cly[0])-1), '1', self.al45y_1, f"{self.k_angle}\n"
        ])
        self.crosslink_bonded['angles'].append([
            str(int(clx[0])-1), clx[0], cly[0], '1', self.al45y_2, f"{self.k_angle}\n"
        ])