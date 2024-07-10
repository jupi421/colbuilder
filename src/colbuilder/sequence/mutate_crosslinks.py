# src/colbuilder/sequence/mutate_crosslinks.py
import os
import logging
from modeller import *
from modeller.scripts import complete_pdb

logger = logging.getLogger(__name__)

def rename_residue_in_pdb(pdb_file: str, chain_id: str, original_resnum: int, new_resname: str):
    with open(pdb_file, "r") as file:
        lines = file.readlines()
    with open(pdb_file, "w") as file:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_num = int(line[22:26].strip())
                chain = line[21].strip()
                if res_num == original_resnum and chain == chain_id:
                    line = line[:17] + new_resname + line[20:]
            file.write(line)

def parse_crosslink_info(crosslink_row):
    if crosslink_row is None:
        return []
    parsed_crosslinks = []
    for i in range(1, 4):  # R1/P1, R2/P2, R3/P3
        patch_type = crosslink_row[f'R{i}']
        if patch_type != 'NONE':
            position = crosslink_row[f'P{i}']
            chain_id, resnum_str = position.split(".")
            resnum = int(resnum_str)
            parsed_crosslinks.append((patch_type, chain_id, resnum))
    return parsed_crosslinks

def apply_crosslinks(input_pdb: str, output_pdb: str, n_crosslink, c_crosslink, restyp_lib: str, top_heav_lib: str, par_mod_lib: str):
    if n_crosslink is None and c_crosslink is None:
        logger.warning("No crosslinks provided. Skipping crosslink application.")
        return

    try:
        env = Environ(rand_seed=-8123, restyp_lib_file=restyp_lib, copy=None)
        env.io.atom_files_directory = ["."]
        env.io.hetatm = True
        env.libs.topology.read(top_heav_lib)
        env.libs.parameters.read(par_mod_lib)

        def patches(mdl):
            mdl.rename_segments(segment_ids=["A", "B", "C"], renumber_residues=[1, 1, 1])
          
            all_patches = parse_crosslink_info(n_crosslink) + parse_crosslink_info(c_crosslink)
              
            for patch_type, chain_id, resnum in all_patches:
                target_residue = f"{resnum}:{chain_id}"
                logger.info(f"Applying {patch_type} patch to target residue: {resnum} in chain {chain_id}")
                mdl.patch(residue_type=patch_type, residues=(mdl.residues[target_residue],))

        mdl = complete_pdb(env, input_pdb, special_patches=patches)
        mdl.write(file=output_pdb)

        # Post-process the output PDB to rename the residues
        all_patches = parse_crosslink_info(n_crosslink) + parse_crosslink_info(c_crosslink)
        for patch_type, chain_id, resnum in all_patches:
            rename_residue_in_pdb(output_pdb, chain_id, resnum, patch_type)

        logger.info(f"Crosslinked PDB saved as: {output_pdb}")
    except Exception as e:
        logger.error(f"An error occurred while applying crosslinks: {str(e)}")
        raise
