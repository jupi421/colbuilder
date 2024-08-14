import chimera
from chimera import runCommand
import os
import random
import math
import sys
import json

residue_rotations = {}

def calculate_center_of_mass(residue):
    total_mass = 0.0
    weighted_sum = chimera.Vector(0.0, 0.0, 0.0)
    for atom in residue.atoms:
        mass = atom.element.mass
        total_mass += mass
        weighted_sum += chimera.Vector(atom.coord().x * mass, atom.coord().y * mass, atom.coord().z * mass)
    return weighted_sum / total_mass

def vector_length(v):
    return math.sqrt(v.x**2 + v.y**2 + v.z**2)

def normalize_vector(v):
    length = vector_length(v)
    return chimera.Vector(v.x/length, v.y/length, v.z/length)

def rotate_residue(residue, axis, angle, center):
    rot = chimera.Xform.rotation(axis, angle)
    for atom in residue.atoms:
        vec = atom.coord() - center
        rotated_vec = rot.apply(vec)
        new_coord = rotated_vec + center
        atom.setCoord(new_coord)
    
    residue_key = (residue.id.chainId, residue.id.position)
    if residue_key not in residue_rotations:
        residue_rotations[residue_key] = []
    residue_rotations[residue_key].append((axis, angle))

def apply_tracked_rotation(residue):
    residue_key = (residue.id.chainId, residue.id.position)
    if residue_key in residue_rotations:
        center = calculate_center_of_mass(residue)
        for axis, angle in residue_rotations[residue_key]:
            rot = chimera.Xform.rotation(axis, angle)
            for atom in residue.atoms:
                vec = atom.coord() - center
                rotated_vec = rot.apply(vec)
                new_coord = rotated_vec + center
                atom.setCoord(new_coord)

def get_rotatable_bonds(residue):
    rotatable_bonds = []
    for atom in residue.atoms:
        if atom.name.startswith('C') and len(atom.bonds) == 1:
            bond = atom.bonds[0]
            if bond.otherAtom(atom).name.startswith('C'):
                rotatable_bonds.append(bond)
    return rotatable_bonds

def rotate_bond(bond, angle):
    atom1, atom2 = bond.atoms
    axis = atom2.coord() - atom1.coord()
    axis = normalize_vector(axis)
    rot = chimera.Xform.rotation(axis, angle)
    movable_atoms = set(atom2.residue.atoms)
    for a in atom2.bondsMap:
        movable_atoms.update(a.residue.atoms)
    movable_atoms.remove(atom1)
    for atom in movable_atoms:
        old_pos = atom.coord()
        new_pos = rot.apply(old_pos - atom1.coord()) + atom1.coord()
        atom.setCoord(new_pos)

def optimize_crosslink_residues(model1, model2, atom1, atom2, is_first_pair, max_distance, max_iterations, initial_temperature, cooling_rate, max_no_improvement):
    if is_first_pair:
        initial_distance = chimera.distance(atom1.coord(), model2.openState.xform.apply(atom2.coord()))
    else:
        initial_distance = chimera.distance(model2.openState.xform.apply(atom1.coord()), atom2.coord())

    best_distance = initial_distance
    best_coords1 = [atom.coord() for atom in atom1.residue.atoms]
    best_coords2 = [atom.coord() for atom in atom2.residue.atoms]

    rotatable_bonds1 = get_rotatable_bonds(atom1.residue)
    rotatable_bonds2 = get_rotatable_bonds(atom2.residue)

    temperature = initial_temperature
    no_improvement_count = 0

    for iteration in range(max_iterations):
        if random.choice([True, False]):
            if rotatable_bonds1:
                bond = random.choice(rotatable_bonds1)
                rotate_bond(bond, random.uniform(-30, 30))
        else:
            if rotatable_bonds2:
                bond = random.choice(rotatable_bonds2)
                rotate_bond(bond, random.uniform(-30, 30))

        if is_first_pair:
            distance = chimera.distance(atom1.coord(), model2.openState.xform.apply(atom2.coord()))
        else:
            distance = chimera.distance(model2.openState.xform.apply(atom1.coord()), atom2.coord())
        
        if distance < best_distance:
            best_distance = distance
            best_coords1 = [atom.coord() for atom in atom1.residue.atoms]
            best_coords2 = [atom.coord() for atom in atom2.residue.atoms]
            no_improvement_count = 0
        else:
            no_improvement_count += 1
            if random.random() < math.exp((best_distance - distance) / temperature):
                best_distance = distance
                best_coords1 = [atom.coord() for atom in atom1.residue.atoms]
                best_coords2 = [atom.coord() for atom in atom2.residue.atoms]
            else:
                for atom, coord in zip(atom1.residue.atoms, best_coords1):
                    atom.setCoord(coord)
                for atom, coord in zip(atom2.residue.atoms, best_coords2):
                    atom.setCoord(coord)

        temperature *= cooling_rate

        if no_improvement_count >= max_no_improvement:
            for residue, bonds in [(atom1.residue, rotatable_bonds1), (atom2.residue, rotatable_bonds2)]:
                if bonds:
                    bond = random.choice(bonds)
                    rotate_bond(bond, random.uniform(-60, 60))
            no_improvement_count = 0
            temperature = min(temperature * 2, initial_temperature)

    for atom, coord in zip(atom1.residue.atoms, best_coords1):
        atom.setCoord(coord)
    for atom, coord in zip(atom2.residue.atoms, best_coords2):
        atom.setCoord(coord)

    return best_distance <= max_distance

def create_collagen_dimer(pdb_file, output_file, crosslink_info, optimization_params):
    original = chimera.openModels.open(pdb_file)[0]
    copy = chimera.openModels.open(pdb_file)[0]
    
    residues_original = {
        crosslink['residue1_type']: sorted(
            [r for r in original.residues if r.type == crosslink['residue1_type']],
            key=lambda r: int(r.id.position)
        )
        for crosslink in crosslink_info
    }
    residues_original.update({
        crosslink['residue2_type']: sorted(
            [r for r in original.residues if r.type == crosslink['residue2_type']],
            key=lambda r: int(r.id.position)
        )
        for crosslink in crosslink_info
    })
    residues_copy = {
        crosslink['residue1_type']: sorted(
            [r for r in copy.residues if r.type == crosslink['residue1_type']],
            key=lambda r: int(r.id.position)
        )
        for crosslink in crosslink_info
    }
    residues_copy.update({
        crosslink['residue2_type']: sorted(
            [r for r in copy.residues if r.type == crosslink['residue2_type']],
            key=lambda r: int(r.id.position)
        )
        for crosslink in crosslink_info
    })
    
    pairs = []
    for i, crosslink in enumerate(crosslink_info):
        if i == 0:
            residue1 = next((r for r in residues_original[crosslink['residue1_type']] if r.id.position == int(crosslink['residue1_position'])), None)
            residue2 = next((r for r in residues_copy[crosslink['residue2_type']] if r.id.position == int(crosslink['residue2_position'])), None)
        else:
            residue1 = next((r for r in residues_copy[crosslink['residue1_type']] if r.id.position == int(crosslink['residue1_position'])), None)
            residue2 = next((r for r in residues_original[crosslink['residue2_type']] if r.id.position == int(crosslink['residue2_position'])), None)
        
        if residue1 is None or residue2 is None:
            print("Error: Could not find specified residues for pair {}".format(i+1))
            return
        
        pairs.append((residue1, residue2, crosslink['atom1'], crosslink['atom2']))
    
    if len(pairs) != 2:
        print("Error: Expected 2 crosslink pairs, but found {}".format(len(pairs)))
        return

    residue1, residue2, atom1_name, atom2_name = pairs[0]
    atom1 = residue1.findAtom(atom1_name)
    atom2 = residue2.findAtom(atom2_name)

    target_distance = optimization_params['target_distance']
    current_vector = atom2.coord() - atom1.coord()
    current_distance = vector_length(current_vector)
    normalized_vector = normalize_vector(current_vector)
    translation_vector = normalized_vector * (current_distance - target_distance)
    
    for atom in copy.atoms:
        new_coord = atom.coord() - translation_vector
        atom.setCoord(new_coord)

    residue3, residue4, atom3_name, atom4_name = pairs[1]
    atom3 = residue3.findAtom(atom3_name)
    atom4 = residue4.findAtom(atom4_name)
    
    second_pair_distance = chimera.distance(atom3.coord(), atom4.coord())

    if second_pair_distance > target_distance:
        current_vector2 = atom4.coord() - atom3.coord()
        current_distance2 = vector_length(current_vector2)
        normalized_vector2 = normalize_vector(current_vector2)
        translation_vector2 = normalized_vector2 * (current_distance2 - target_distance) * 0.5
        
        for atom in copy.atoms:
            new_coord = atom.coord() - translation_vector2
            atom.setCoord(new_coord)

    for i, (res1, res2, atom1_name, atom2_name) in enumerate(pairs):
        atom1 = res1.findAtom(atom1_name)
        atom2 = res2.findAtom(atom2_name)
        
        optimize_crosslink_residues(original, copy, atom1, atom2, is_first_pair=(i==0), 
                                    max_distance=optimization_params['max_distance'], 
                                    max_iterations=optimization_params['max_iterations'],
                                    initial_temperature=optimization_params['initial_temperature'],
                                    cooling_rate=optimization_params['cooling_rate'],
                                    max_no_improvement=optimization_params['max_no_improvement'])
        
        apply_tracked_rotation(res1)
        apply_tracked_rotation(res2)
        
        # Update the corresponding residues in the original model
        original_res1 = next(r for r in original.residues if r.id == res1.id)
        original_res2 = next(r for r in original.residues if r.id == res2.id)
        
        if res1 in residues_original[res1.type]:
            # res1 is from the original model, just update coordinates
            for orig_atom, rotated_atom in zip(original_res1.atoms, res1.atoms):
                orig_atom.setCoord(rotated_atom.coord())
        else:
            # res1 is from the copy, apply the same rotations to the original coordinates
            center = calculate_center_of_mass(original_res1)
            for axis, angle in residue_rotations.get((res1.id.chainId, res1.id.position), []):
                rot = chimera.Xform.rotation(axis, angle)
                for atom in original_res1.atoms:
                    vec = atom.coord() - center
                    rotated_vec = rot.apply(vec)
                    new_coord = rotated_vec + center
                    atom.setCoord(new_coord)

        if res2 in residues_original[res2.type]:
            # res2 is from the original model, just update coordinates
            for orig_atom, rotated_atom in zip(original_res2.atoms, res2.atoms):
                orig_atom.setCoord(rotated_atom.coord())
        else:
            # res2 is from the copy, apply the same rotations to the original coordinates
            center = calculate_center_of_mass(original_res2)
            for axis, angle in residue_rotations.get((res2.id.chainId, res2.id.position), []):
                rot = chimera.Xform.rotation(axis, angle)
                for atom in original_res2.atoms:
                    vec = atom.coord() - center
                    rotated_vec = rot.apply(vec)
                    new_coord = rotated_vec + center
                    atom.setCoord(new_coord)

    runCommand("write format pdb 0 {}".format(output_file))

if __name__ == "__main__":
    input_pdb = os.environ.get('INPUT_PDB')
    output_pdb = os.environ.get('OUTPUT_PDB')
    crosslink_info_json = os.environ.get('CROSSLINK_INFO')
    optimization_params_json = os.environ.get('OPTIMIZATION_PARAMS')

    if input_pdb and output_pdb and crosslink_info_json and optimization_params_json:
        try:
            crosslink_info = json.loads(crosslink_info_json)
            optimization_params = json.loads(optimization_params_json)
            create_collagen_dimer(input_pdb, output_pdb, crosslink_info, optimization_params)
        except Exception as e:
            sys.stderr.write("An error occurred in the Chimera script: {}\n".format(str(e)))
            raise
    else:
        sys.stderr.write("Error: Required environment variables not set\n")