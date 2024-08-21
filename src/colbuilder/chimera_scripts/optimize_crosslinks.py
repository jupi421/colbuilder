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

    print("Starting optimization for atoms {} and {}".format(atom1, atom2))
    print("Initial distance: {:.3f}".format(initial_distance))

    temperature = params['initial_temperature']
    for iteration in range(params['max_iterations']):
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
            print("Iteration {}: New best distance: {:.3f}".format(iteration, best_distance))
        else:
            for atom, coord in zip(atom1.residue.atoms, best_coords1):
                atom.setCoord(coord)
            for atom, coord in zip(atom2.residue.atoms, best_coords2):
                atom.setCoord(coord)

        temperature *= params['cooling_rate']
        
        if iteration % 100 == 0:
            print("Iteration {}: Current distance: {:.3f}, Best distance: {:.3f}, Temperature: {:.2f}".format(
                iteration, current_distance, best_distance, temperature))
        
        if best_distance <= params['target_distance']:
            print("Target distance reached at iteration {}".format(iteration))
            break

    print("Optimization completed. Final distance: {:.3f}".format(best_distance))
    
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

def create_optimized_collagen(input_pdbs, output_pdb, crosslink_info, optimization_params):
    models = []
    for pdb in input_pdbs:
        models.extend(chimera.openModels.open(pdb))
    
    print("Loaded {} models for optimization.".format(len(models)))

    original_model = models[0]  
    transformations_to_apply = {}  

    for i, crosslink in enumerate(crosslink_info):
        print("\nProcessing crosslink {} of {}".format(i+1, len(crosslink_info)))
        
        if crosslink['residue3_type'] != "NONE":
            # Trivalent crosslink
            residue1_info = {'type': str(crosslink['residue1_type']), 'position': int(crosslink['residue1_position'])}
            residue2_info = {'type': str(crosslink['residue2_type']), 'position': int(crosslink['residue2_position'])}
            residue3_info = {'type': str(crosslink['residue3_type']), 'position': int(crosslink['residue3_position'])}
            
            closest_pair1, initial_distance1 = find_closest_pair(models, residue1_info, residue3_info, str(crosslink['atom1']), str(crosslink['atom31']))
            closest_pair2, initial_distance2 = find_closest_pair(models, residue2_info, residue3_info, str(crosslink['atom2']), str(crosslink['atom32']))
            
            if not (closest_pair1 and closest_pair2):
                print("Warning: Couldn't find suitable pairs for trivalent crosslink. Skipping.")
                continue

            atom1, atom31 = closest_pair1
            atom2, atom32 = closest_pair2

            print("Proceeding with optimization for trivalent crosslink:")
            print("  Pair 1: {} {} - {} {}".format(atom1.residue.type, atom1.residue.id.position, atom31.residue.type, atom31.residue.id.position))
            print("  Initial distance 1: {:.3f}".format(initial_distance1))
            print("  Pair 2: {} {} - {} {}".format(atom2.residue.type, atom2.residue.id.position, atom32.residue.type, atom32.residue.id.position))
            print("  Initial distance 2: {:.3f}".format(initial_distance2))

            if initial_distance1 > 100 or initial_distance2 > 100:
                print("Warning: Initial distances are unusually large. Skipping this crosslink.")
                continue

            success1, transformations1 = optimize_crosslink(atom1, atom31, optimization_params)
            success2, transformations2 = optimize_crosslink(atom2, atom32, optimization_params)

            if success1 and success2:
                transformations_to_apply[(atom1.residue.type, atom1.residue.id.position)] = transformations1[0]
                transformations_to_apply[(atom31.residue.type, atom31.residue.id.position)] = transformations1[1]
                transformations_to_apply[(atom2.residue.type, atom2.residue.id.position)] = transformations2[0]
                # Note: atom32 is the same residue as atom31, so we don't add it again
            else:
                print("Warning: Optimization failed for this trivalent crosslink.")

        else:
            # Divalent crosslink
            residue1_info = {'type': str(crosslink['residue1_type']), 'position': int(crosslink['residue1_position'])}
            residue2_info = {'type': str(crosslink['residue2_type']), 'position': int(crosslink['residue2_position'])}
            
            closest_pair, initial_distance = find_closest_pair(models, residue1_info, residue2_info, str(crosslink['atom1']), str(crosslink['atom2']))
            
            if not closest_pair:
                print("Warning: Couldn't find a suitable pair for divalent crosslink. Skipping.")
                continue

            atom1, atom2 = closest_pair
            print("Proceeding with optimization for divalent crosslink:")
            print("  {} {} - {} {}".format(atom1.residue.type, atom1.residue.id.position, atom2.residue.type, atom2.residue.id.position))
            print("  Initial distance: {:.3f}".format(initial_distance))

            if initial_distance > 100:
                print("Warning: Initial distance is unusually large. Skipping this pair.")
                continue

            success, transformations = optimize_crosslink(atom1, atom2, optimization_params)
            if success:
                transformations_to_apply[(atom1.residue.type, atom1.residue.id.position)] = transformations[0]
                transformations_to_apply[(atom2.residue.type, atom2.residue.id.position)] = transformations[1]
            else:
                print("Warning: Optimization failed for this divalent crosslink.")

    for residue in original_model.residues:
        residue_key = (residue.type, residue.id.position)
        if residue_key in transformations_to_apply:
            apply_transformation(residue, transformations_to_apply[residue_key])
            print("Applied transformation to residue {} {}".format(residue.type, residue.id.position))

    runCommand("write format pdb #{} {}".format(original_model.id, output_pdb))

    print("Optimization completed. Output saved to: {}".format(output_pdb))

    chimera.openModels.close(models[1:])

def main():
    try:
        input_pdbs = os.environ.get('INPUT_PDB', '').split()
        output_pdb = os.environ.get('OUTPUT_PDB', '')
        crosslink_info_json = os.environ.get('CROSSLINK_INFO', '')
        optimization_params_json = os.environ.get('OPTIMIZATION_PARAMS', '')

        if all([input_pdbs, output_pdb, crosslink_info_json, optimization_params_json]):
            print("Starting collagen crosslink optimization")
            print("Input PDBs: {}".format(input_pdbs))
            print("Output PDB: {}".format(output_pdb))
            crosslink_info = json.loads(crosslink_info_json)
            optimization_params = json.loads(optimization_params_json)
            print("Crosslink info: {}".format(json.dumps(crosslink_info, indent=2)))
            print("Optimization parameters: {}".format(json.dumps(optimization_params, indent=2)))
            create_optimized_collagen(input_pdbs, output_pdb, crosslink_info, optimization_params)
            print("Optimization completed successfully")
        else:
            raise ValueError("Required environment variables not set")
    except Exception as e:
        sys.stderr.write("An error occurred in the Chimera script: {}\n".format(str(e)))
        chimera.closeSession()
        sys.exit(1)
    finally:
        chimera.closeSession()
        sys.exit(0)

if __name__ == "__main__":
    main()