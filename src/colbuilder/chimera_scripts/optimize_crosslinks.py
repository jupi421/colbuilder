import chimera
from chimera import runCommand
import os
import random
import math
import sys
import json

try:
    chimera.dialogs.find("reply", create=False).Close()
except:
    pass

log_file = os.path.join(os.path.expanduser("~"), "chimera_script.log")
def log(message):
    with open(log_file, "a") as f:
        f.write(message + "\n")

log("Chimera script started")

def vector_length(v):
    return math.sqrt(v.x**2 + v.y**2 + v.z**2)

def normalize_vector(v):
    length = vector_length(v)
    return chimera.Vector(v.x/length, v.y/length, v.z/length)

def rotate_residue(residue, axis, angle):
    center = chimera.Vector(0, 0, 0)
    for atom in residue.atoms:
        center += chimera.Vector(atom.coord().x, atom.coord().y, atom.coord().z)
    center /= len(residue.atoms)
    
    rotation = chimera.Xform.rotation(axis, angle)
    for atom in residue.atoms:
        vec = chimera.Vector(atom.coord().x, atom.coord().y, atom.coord().z) - center
        rotated_vec = rotation.apply(vec)
        new_coord = rotated_vec + center
        atom.setCoord(chimera.Coord(new_coord.x, new_coord.y, new_coord.z))

def find_closest_pair(models, residue1_info, residue2_info, atom1_name, atom2_name):
    closest_pair = None
    min_distance = float('inf')
    
    print("Searching for closest pair:")
    print("  Residue 1: {} {}, Atom: {}".format(residue1_info['type'], residue1_info['position'], atom1_name))
    print("  Residue 2: {} {}, Atom: {}".format(residue2_info['type'], residue2_info['position'], atom2_name))
    
    for m1 in models:
        residues1 = [r for r in m1.residues if str(r.type) == str(residue1_info['type']) and r.id.position == int(residue1_info['position'])]
        print("\nModel #{}: Found {} matching residues for Residue 1 ({} {}).".format(m1.id, len(residues1), residue1_info['type'], residue1_info['position']))
        
        for r1 in residues1:
            atom1 = r1.findAtom(str(atom1_name))
            if not atom1:
                print("  Warning: Atom {} not found in residue {} {} of model #{}".format(atom1_name, r1.type, r1.id.position, m1.id))
                continue
            
            print("  Residue 1: {} {} (model #{}) with atom {} at coordinates {}".format(r1.type, r1.id.position, m1.id, atom1_name, atom1.coord()))
            
            for m2 in models:
                residues2 = [r for r in m2.residues if str(r.type) == str(residue2_info['type']) and r.id.position == int(residue2_info['position'])]
                print("  Model #{}: Found {} matching residues for Residue 2 ({} {}).".format(m2.id, len(residues2), residue2_info['type'], residue2_info['position']))
                
                for r2 in residues2:
                    atom2 = r2.findAtom(str(atom2_name))
                    if not atom2:
                        print("    Warning: Atom {} not found in residue {} {} of model #{}".format(atom2_name, r2.type, r2.id.position, m2.id))
                        continue
                    
                    distance = chimera.distance(atom1.coord(), atom2.coord())
                    print("    Residue 2: {} {} (model #{}) with atom {} at coordinates {}".format(r2.type, r2.id.position, m2.id, atom2_name, atom2.coord()))
                    print("    Distance between Residue 1 ({} {} model #{}) and Residue 2 ({} {} model #{}): {:.3f}".format(
                        r1.type, r1.id.position, m1.id,
                        r2.type, r2.id.position, m2.id,
                        distance))
                    
                    if distance < min_distance:
                        min_distance = distance
                        closest_pair = (atom1, atom2, m1.id, m2.id)
                        print("    New closest pair found!")

    if closest_pair:
        print("\nClosest pair found:")
        print("  Atom 1: {} {} (model #{})".format(closest_pair[0].residue.type, closest_pair[0].residue.id.position, closest_pair[2]))
        print("  Atom 2: {} {} (model #{})".format(closest_pair[1].residue.type, closest_pair[1].residue.id.position, closest_pair[3]))
        print("  Distance: {:.3f}".format(min_distance))
        runCommand("show #{}:{}@{} #{}:{}@{}".format(
            closest_pair[2], closest_pair[0].residue.id.position, closest_pair[0].name,
            closest_pair[3], closest_pair[1].residue.id.position, closest_pair[1].name))
        #runCommand("write format pdb #{} closest_pair.pdb".format(closest_pair[2]))
    else:
        print("\nNo suitable pair found.")

    return closest_pair[:2], min_distance


def optimize_crosslink(atom1, atom2, params):
    initial_distance = chimera.distance(atom1.coord(), atom2.coord())
    best_distance = initial_distance
    best_coords1 = [atom.coord() for atom in atom1.residue.atoms]
    best_coords2 = [atom.coord() for atom in atom2.residue.atoms]
    transformations_applied1 = []
    transformations_applied2 = []

    log("Starting optimization for atoms {} and {}".format(atom1, atom2))
    log("Initial distance: {:.3f}".format(initial_distance))

    temperature = params['initial_temperature']
    for iteration in xrange(params['max_iterations']):
        axis1 = normalize_vector(chimera.Vector(*(atom1.residue.atoms[1].coord() - atom1.residue.atoms[0].coord())))
        axis2 = normalize_vector(chimera.Vector(*(atom2.residue.atoms[1].coord() - atom2.residue.atoms[0].coord())))
        angle1 = random.uniform(-5, 5)
        angle2 = random.uniform(-5, 5)
        rotation1 = chimera.Xform.rotation(axis1, angle1)
        rotation2 = chimera.Xform.rotation(axis2, angle2)
        
        rotate_residue(atom1.residue, axis1, angle1)
        rotate_residue(atom2.residue, axis2, angle2)

        current_distance = chimera.distance(atom1.coord(), atom2.coord())
        
        if current_distance < best_distance or random.random() < math.exp((best_distance - current_distance) / temperature):
            best_distance = current_distance
            best_coords1 = [atom.coord() for atom in atom1.residue.atoms]
            best_coords2 = [atom.coord() for atom in atom2.residue.atoms]
            transformations_applied1.append((axis1, angle1))
            transformations_applied2.append((axis2, angle2))
            log("Iteration {}: New best distance: {:.3f}".format(iteration, best_distance))
        else:
            for atom, coord in zip(atom1.residue.atoms, best_coords1):
                atom.setCoord(coord)
            for atom, coord in zip(atom2.residue.atoms, best_coords2):
                atom.setCoord(coord)

        temperature *= params['cooling_rate']
        
        if iteration % 100 == 0:
            log("Iteration {}: Current distance: {:.3f}, Best distance: {:.3f}, Temperature: {:.2f}".format(
                iteration, current_distance, best_distance, temperature))
        
        if best_distance <= params['target_distance']:
            log("Target distance reached at iteration {}".format(iteration))
            break

    log("Optimization completed. Final distance: {:.3f}".format(best_distance))
    
    return best_distance <= params['target_distance'], (transformations_applied1, transformations_applied2)

def apply_transformation(residue, transformations):
    center = chimera.Vector(0, 0, 0)
    for atom in residue.atoms:
        center += chimera.Vector(atom.coord().x, atom.coord().y, atom.coord().z)
    center /= len(residue.atoms)
    
    for axis, angle in transformations:
        rotation = chimera.Xform.rotation(axis, angle)
        for atom in residue.atoms:
            vec = atom.coord() - center
            new_coord = rotation.apply(vec) + center
            atom.setCoord(chimera.Coord(new_coord.x, new_coord.y, new_coord.z))

    log("Applied {} transformations to residue {} {}".format(len(transformations), residue.type, residue.id.position))

def create_optimized_collagen(input_pdbs, output_pdb, crosslink_info, optimization_params):
    models = []
    for pdb in input_pdbs:
        models.extend(chimera.openModels.open(pdb))
    
    log("Loaded {} models for optimization.".format(len(models)))

    original_model = models[0]  
    transformations_to_apply = {}  

    for i, crosslink in enumerate(crosslink_info):
        log("\nProcessing crosslink {} of {}".format(i+1, len(crosslink_info)))
        
        residue1_info = {'type': str(crosslink['residue1_type']), 'position': int(crosslink['residue1_position'])}
        residue2_info = {'type': str(crosslink['residue2_type']), 'position': int(crosslink['residue2_position'])}
        
        closest_pair, initial_distance = find_closest_pair(models, residue1_info, residue2_info, str(crosslink['atom1']), str(crosslink['atom2']))
        
        if not closest_pair:
            log("Warning: Couldn't find a suitable pair for crosslink. Skipping.")
            continue

        atom1, atom2 = closest_pair
        log("Proceeding with optimization for:")
        log("  Atom 1: {} {} (model #{})".format(atom1.residue.type, atom1.residue.id.position, atom1.molecule.id))
        log("  Atom 2: {} {} (model #{})".format(atom2.residue.type, atom2.residue.id.position, atom2.molecule.id))
        log("  Initial distance: {:.3f}".format(initial_distance))

        if initial_distance > 100:  
            log("Warning: Initial distance is unusually large. Skipping this pair.")
            continue

        success, transformations = optimize_crosslink(atom1, atom2, optimization_params)
        if not success:
            log("Warning: Optimization failed for this crosslink.")
        else:
            residue1_key = (atom1.residue.type, atom1.residue.id.position)
            residue2_key = (atom2.residue.type, atom2.residue.id.position)
            transformations_to_apply[residue1_key] = transformations[0]
            transformations_to_apply[residue2_key] = transformations[1]

    for residue in original_model.residues:
        residue_key = (residue.type, residue.id.position)
        if residue_key in transformations_to_apply:
            apply_transformation(residue, transformations_to_apply[residue_key])

    log("\nSaving the optimized original model...")
    runCommand("write format pdb 0 {}".format(output_pdb))

    log("Closing additional models...")
    for m in models[1:]:
        chimera.openModels.close(m)

    log("Optimization completed. Output saved to: {}".format(output_pdb))

if __name__ == "__main__":
    input_pdbs = os.environ.get('INPUT_PDB').split()
    output_pdb = os.environ.get('OUTPUT_PDB')
    crosslink_info_json = os.environ.get('CROSSLINK_INFO')
    optimization_params_json = os.environ.get('OPTIMIZATION_PARAMS')

    if all([input_pdbs, output_pdb, crosslink_info_json, optimization_params_json]):
        try:
            print("Starting collagen crosslink optimization")
            print("Input PDBs: {}".format(input_pdbs))
            print("Output PDB: {}".format(output_pdb))
            crosslink_info = json.loads(crosslink_info_json)
            optimization_params = json.loads(optimization_params_json)
            print("Crosslink info: {}".format(json.dumps(crosslink_info, indent=2)))
            print("Optimization parameters: {}".format(json.dumps(optimization_params, indent=2)))
            create_optimized_collagen(input_pdbs, output_pdb, crosslink_info, optimization_params)
        except Exception as e:
            sys.stderr.write("An error occurred in the Chimera script: {}\n".format(str(e)))
            raise
    else:
        sys.stderr.write("Error: Required environment variables not set\n")
