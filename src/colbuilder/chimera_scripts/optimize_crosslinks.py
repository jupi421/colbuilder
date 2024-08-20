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

def rotate_ca_cb_bond(residue, angle):
    ca_atom = residue.findAtom('CA')
    cb_atom = residue.findAtom('CB')
    if not ca_atom or not cb_atom:
        return None

    axis = cb_atom.coord() - ca_atom.coord()
    axis = normalize_vector(axis)
    rot = chimera.Xform.rotation(axis, angle)

    for atom in residue.atoms:
        if atom.name not in ['CA', 'CB']:
            vec = atom.coord() - ca_atom.coord()
            rotated_vec = rot.apply(vec)
            new_coord = rotated_vec + ca_atom.coord()
            atom.setCoord(new_coord)
    
    residue_key = (residue.id.chainId, residue.id.position)
    if residue_key not in residue_rotations:
        residue_rotations[residue_key] = []
    residue_rotations[residue_key].append((axis, angle))

def translate_residue(residue, translation_vector):
    for atom in residue.atoms:
        new_coord = atom.coord() + translation_vector
        atom.setCoord(new_coord)
    return ('translation', translation_vector)

def optimize_residues_by_rotation(atom1, atom2, max_distance, max_iterations, temperature, step_size_reduction, max_no_improvement, translation_step_size):
    initial_distance = chimera.distance(atom1.coord(), atom2.coord())
    # print("Initial distance: {:.2f} Angstroms".format(initial_distance))

    best_distance = initial_distance
    best_coords1 = [atom.coord() for atom in atom1.residue.atoms]
    best_coords2 = [atom.coord() for atom in atom2.residue.atoms]

    no_improvement_count = 0
    step_size = 5  
    min_step_size = 1  
    translation_attempted = False

    for iteration in range(max_iterations):
        # Rotate CA-CB bonds with small, controlled angles
        rotation1 = rotate_ca_cb_bond(atom1.residue, random.uniform(-step_size, step_size))
        rotation2 = rotate_ca_cb_bond(atom2.residue, random.uniform(-step_size, step_size))

        current_distance = chimera.distance(atom1.coord(), atom2.coord())
        # print("Current distance: {:.2f} Angstroms".format(current_distance))

        delta_distance = current_distance - best_distance

        if delta_distance < 0 or random.random() < math.exp(-delta_distance / temperature):
            # print("Transformation accepted. Delta distance: {:.2f}".format(delta_distance))
            best_distance = current_distance
            best_coords1 = [atom.coord() for atom in atom1.residue.atoms]
            best_coords2 = [atom.coord() for atom in atom2.residue.atoms]
            no_improvement_count = 0
            step_size = max(min_step_size, step_size * step_size_reduction)
            translation_attempted = False
        else:
            # print("Transformation rejected. Reverting rotation.")
            for atom, coord in zip(atom1.residue.atoms, best_coords1):
                atom.setCoord(coord)
            for atom, coord in zip(atom2.residue.atoms, best_coords2):
                atom.setCoord(coord)
            no_improvement_count += 1

        if no_improvement_count > max_no_improvement:
            if not translation_attempted:
                # print("No improvement. Attempting translation.")
                direction = normalize_vector(atom2.coord() - atom1.coord())
                translation_vector = direction * translation_step_size
                translate_residue(atom1.residue, translation_vector)
                translate_residue(atom2.residue, -translation_vector)
                translation_attempted = True
                no_improvement_count = 0
            else:
                # print("No improvement. Increasing step size.")
                step_size *= 1.5
                translation_attempted = False

        if best_distance <= max_distance:
            # print("Acceptable distance reached. Stopping optimization.")
            break

        temperature *= 0.99

    for atom, coord in zip(atom1.residue.atoms, best_coords1):
        atom.setCoord(coord)
    for atom, coord in zip(atom2.residue.atoms, best_coords2):
        atom.setCoord(coord)

    # print("Final best distance: {:.2f} Angstroms".format(best_distance))
    return best_distance <= max_distance

def apply_transformations_to_original(original_residue, transformed_residue, residue_rotations):
    center = calculate_center_of_mass(original_residue)
    for axis, angle in residue_rotations.get((transformed_residue.id.chainId, transformed_residue.id.position), []):
        rot = chimera.Xform.rotation(axis, angle)
        for atom in original_residue.atoms:
            vec = atom.coord() - center
            rotated_vec = rot.apply(vec)
            new_coord = rotated_vec + center
            atom.setCoord(new_coord)

def save_dimer(original, copy, dimer_output_file):
    with open(dimer_output_file, 'w') as f:
        runCommand("write format pdb 0 output_original.pdb")
        with open("output_original.pdb", 'r') as original_file:
            for line in original_file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    f.write(line)

        runCommand("write format pdb 1 output_copy.pdb")
        with open("output_copy.pdb", 'r') as copy_file:
            for line in copy_file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    f.write(line)

def create_collagen_dimer(pdb_file, output_file, dimer_output_file, crosslink_info, optimization_params):
    # print("Starting collagen dimer creation")
    # print("Input PDB file: {}".format(pdb_file))
    # print("Output PDB file: {}".format(output_file))
    # print("Dimer output PDB file: {}".format(dimer_output_file))
    # print("Crosslink info: {}".format(json.dumps(crosslink_info, indent=2)))
    # print("Optimization parameters: {}".format(json.dumps(optimization_params, indent=2)))

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

    current_vector2 = atom4.coord() - atom3.coord()
    current_distance2 = vector_length(current_vector2)
    normalized_vector2 = normalize_vector(current_vector2)
    translation_vector2 = normalized_vector2 * (current_distance2 - target_distance) * 0.5

    original_coords = {atom: atom.coord() for atom in copy.atoms}

    for atom in copy.atoms:
        new_coord = atom.coord() - translation_vector2
        atom.setCoord(new_coord)

    new_distance2 = chimera.distance(atom3.coord(), atom4.coord())
    if new_distance2 > current_distance2:
        for atom, coord in original_coords.items():
            atom.setCoord(coord)
        print("Translation reverted as it increased the distance.")
    else:
        print("Translation applied. New distance: {:.2f} Angstroms".format(new_distance2))

    for i, (res1, res2, atom1_name, atom2_name) in enumerate(pairs):
        atom1 = res1.findAtom(atom1_name)
        atom2 = res2.findAtom(atom2_name)
        # print("\nOptimizing pair {}: {} (Chain {}, Position {}) - {} (Chain {}, Position {})".format(
        #     i+1, res1.type, res1.id.chainId, res1.id.position, res2.type, res2.id.chainId, res2.id.position))

        optimize_residues_by_rotation(atom1, atom2, 
                                      max_distance=optimization_params['max_distance'], 
                                      max_iterations=optimization_params['max_iterations'],
                                      temperature=optimization_params['initial_temperature'], 
                                      step_size_reduction=optimization_params['step_size_reduction'],
                                      max_no_improvement=optimization_params['max_no_improvement'],
                                      translation_step_size=optimization_params['translation_step_size'])

        apply_transformations_to_original(res1, res1, residue_rotations)
        apply_transformations_to_original(res2, res2, residue_rotations)

    runCommand("write format pdb 0 {}".format(output_file)) 

    # # Save both the original and copy models as separate chains in the dimer PDB file
    # save_dimer(original, copy, dimer_output_file)
    # print("\nOptimized model created successfully. Output saved to {} and {}".format(output_file, dimer_output_file))

if __name__ == "__main__":
    input_pdb = os.environ.get('INPUT_PDB')
    output_pdb = os.environ.get('OUTPUT_PDB')
    dimer_output_pdb = "ouput_dimer.pdb"
    crosslink_info_json = os.environ.get('CROSSLINK_INFO')
    optimization_params_json = os.environ.get('OPTIMIZATION_PARAMS')

    if input_pdb and output_pdb and dimer_output_pdb and crosslink_info_json and optimization_params_json:
        try:
            crosslink_info = json.loads(crosslink_info_json)
            optimization_params = json.loads(optimization_params_json)
            create_collagen_dimer(input_pdb, output_pdb, dimer_output_pdb, crosslink_info, optimization_params)
        except Exception as e:
            sys.stderr.write("An error occurred in the Chimera script: {}\n".format(str(e)))
            raise
    else:
        sys.stderr.write("Error: Required environment variables not set\n")