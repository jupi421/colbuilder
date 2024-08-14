"""
This module provides functionality to apply crosslinks to protein structures using MODELLER®.

MODELLER® is a trademark of the Regents of the University of California.
It was developed by Andrej Sali and colleagues at the University of California, San Francisco.
For more information about MODELLER, please visit: https://salilab.org/modeller/

This module is designed to facilitate the use of MODELLER within the Colbuilder project,
but it is not affiliated with or endorsed by the MODELLER developers or the University of California.

When using this module, please ensure you comply with MODELLER's license terms and provide
appropriate attribution in your work.
"""

# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import os
from typing import Optional, List, Tuple
import pandas as pd
from modeller import Environ
from modeller.scripts import complete_pdb

from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.config import ColbuilderConfig

LOG = setup_logger(__name__)

def rename_residue_in_pdb(pdb_file: str, chain_id: str, original_resnum: int, new_resname: str) -> None:
    """
    Rename a residue in a PDB file.

    Args:
        pdb_file (str): Path to PDB file.
        chain_id (str): Chain ID of the residue.
        original_resnum (int): Original residue number.
        new_resname (str): New residue name.

    Raises:
        IOError: If there's an issue reading or writing the PDB file.
    """
    try:
        with open(pdb_file, "r") as file:
            lines = file.readlines()
        with open(pdb_file, "w") as file:
            for line in lines:
                if line.startswith(("ATOM", "HETATM")):
                    res_num = int(line[22:26].strip())
                    chain = line[21].strip()
                    if res_num == original_resnum and chain == chain_id:
                        line = f"{line[:17]}{new_resname}{line[20:]}"
                file.write(line)
    except IOError as e:
        LOG.error(f"Error renaming residue in PDB file: {str(e)}")
        raise

def parse_crosslink_info(crosslink_row: Optional[pd.Series]) -> List[Tuple[str, str, int]]:
    """
    Parse crosslink information from a pandas Series.

    Args:
        crosslink_row (Optional[pd.Series]): A row from the crosslink DataFrame.

    Returns:
        List[Tuple[str, str, int]]: Parsed crosslink information as a list of tuples.
    """
    if crosslink_row is None:
        return []
    
    parsed_crosslinks = []
    for i in range(1, 4):  # R1/P1, R2/P2, R3/P3
        patch_type = crosslink_row.get(f'R{i}', 'NONE')
        if patch_type != 'NONE':
            position = crosslink_row.get(f'P{i}')
            if position and '.' in position:
                resnum_str, chain_id = position.split(".")
                try:
                    resnum = int(resnum_str)
                    parsed_crosslinks.append((patch_type, chain_id, resnum))
                except ValueError:
                    LOG.warning(f"Invalid residue number in position {position} for patch type {patch_type}")
            else:
                LOG.warning(f"Invalid position format: {position} for patch type {patch_type}")
    return parsed_crosslinks

@timeit
def apply_crosslinks(input_pdb: str, output_pdb: str, n_crosslink: Optional[pd.Series], c_crosslink: Optional[pd.Series], cfg: ColbuilderConfig) -> str:
    """
    Apply crosslinks to a PDB file using MODELLER.

    Args:
        input_pdb (str): Path to input PDB file.
        output_pdb (str): Path to output PDB file.
        n_crosslink (Optional[pd.Series]): N-terminal crosslink information.
        c_crosslink (Optional[pd.Series]): C-terminal crosslink information.
        cfg (ColbuilderConfig): Configuration object.

    Returns:
        str: Path to the output PDB file.

    Raises:
        Exception: If an error occurs during crosslink application.
    """
    if n_crosslink is None and c_crosslink is None:
        LOG.warning("No crosslinks provided. Skipping crosslink application.")
        return input_pdb
    
    try:
        LOG.debug(f"N-terminal crosslink: {n_crosslink}")
        LOG.debug(f"C-terminal crosslink: {c_crosslink}")
        
        env = Environ(rand_seed=-8123, restyp_lib_file=str(cfg.RESTYP_LIB_PATH), copy=None)
        env.io.atom_files_directory = ["."]
        env.io.hetatm = True
        env.libs.topology.read(str(cfg.TOP_HEAV_LIB_PATH))
        env.libs.parameters.read(str(cfg.PAR_MOD_LIB_PATH))
        
        def patches(mdl) -> None:
            mdl.rename_segments(segment_ids=["A", "B", "C"], renumber_residues=[1, 1, 1])
            
            all_patches = parse_crosslink_info(n_crosslink) + parse_crosslink_info(c_crosslink)
            LOG.debug(f"All patches to be applied: {all_patches}")

            for patch_type, chain_id, resnum in all_patches:
                target_residue = f"{resnum}:{chain_id}"
                LOG.info(f"   - Applying {patch_type} patch to target residue: {resnum} in chain {chain_id}")
                mdl.patch(residue_type="R" + patch_type, residues=(mdl.residues[target_residue],))

        mdl = complete_pdb(env, input_pdb, special_patches=patches)
        mdl.write(file=output_pdb)

        all_patches = parse_crosslink_info(n_crosslink) + parse_crosslink_info(c_crosslink)
        for patch_type, chain_id, resnum in all_patches:
            rename_residue_in_pdb(output_pdb, chain_id, resnum, patch_type)
        
        LOG.debug(f"Crosslinks applied successfully. Output PDB: {output_pdb}")
        return output_pdb
        
    except Exception as e:
        LOG.error(f"An error occurred while applying crosslinks: {str(e)}")
        raise