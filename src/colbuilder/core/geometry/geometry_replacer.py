"""
Colbuilder Crosslink Replacement Module

This module provides a unified approach to replacing crosslinks with standard amino acids
in a collagen microfibril, supporting both system-based and direct replacement approaches.
"""

import os
import random
import time
import math
import traceback
import subprocess
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, List, Union, Tuple, cast, Set

from ..utils.exceptions import GeometryGenerationError
from ..utils.logger import setup_logger
from ..utils.config import ColbuilderConfig

LOG = setup_logger(__name__)

class CrosslinkReplacer:
    """
    Unified service for replacing crosslinks in collagen structures.
    
    This class consolidates different replacement approaches:
    1. System-based replacement (when working with a System object from geometry generation)
    2. Direct replacement (when working directly with PDB files)
    
    It provides methods for selecting crosslinks to replace, generating
    replacement instructions, running Chimera to perform replacements,
    and handling PDB file operations.
    """
    
    def __init__(self) -> None:
        """
        Initialize the crosslink replacer service.
        
        Sets up the basic state required for the replacer, including
        which line prefixes to process and a container for external data.
        """
        self.is_line = ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
        self.external = {}
        
    async def replace_in_system(self, system: Any, config: ColbuilderConfig) -> Any:
        """
        Replace crosslinks in the provided system according to configuration settings.
        
        This is the main entry point for system-based replacement. It analyzes the system,
        selects crosslinks to replace based on the specified ratio, and performs the
        replacement using Chimera.
        
        Args:
            system: System containing the models to be modified (might be None)
            config: Configuration for replacement including ratio and other settings
            
        Returns:
            Modified system with replaced crosslinks
            
        Raises:
            GeometryGenerationError: If replacement fails for any reason
        """
        try:
            if config.replace_file and os.path.exists(config.replace_file):
                LOG.info(f"Using external replacement file: {config.replace_file}")
                
                with open(config.replace_file, 'r') as f:
                    first_line = f.readline().strip()
                    if first_line.startswith(("ATOM", "CRYST1", "HETATM")):
                        await self.replace_direct(config)
                        return system
                
                if system is None:
                    raise GeometryGenerationError(
                        message="System is None but replacement instructions file was provided. Cannot continue.",
                        error_code="GEO_ERR_004"
                    )
                    
                replace_file = str(config.replace_file)
                await self._run_replacement_with_chimera(system, config, replace_file)
                return system
                
            if not hasattr(system, 'config'):
                system.config = config
            
            z_bounds = self._calculate_fibril_bounds(system, config.fibril_length)
            
            self._select_replacements(system, config.ratio_replace, z_bounds)
            
            replace_file = os.path.join(config.working_directory, "replace.txt")
            self._write_replace_file(system, replace_file)
            
            await self._run_replacement_with_chimera(system, config, replace_file)
            
            replaced_count = system.count_states(state='replace')
            return system
                
        except Exception as e:
            LOG.error(f"Error in system-based replacement: {str(e)}")
            traceback.print_exc()
            raise GeometryGenerationError(
                message=f"Failed to replace crosslinks: {str(e)}",
                original_error=e,
                error_code="GEO_ERR_004",
                context={"config": config.model_dump()}
            )
    
    async def replace_direct(self, config: ColbuilderConfig) -> None:
        """
        Directly replace crosslinks in a PDB file without requiring a System object.
        
        This method provides a standalone replacement approach working directly with PDB files.
        It identifies complete crosslink pairs based on proximity and residue type,
        and ensures both residues in each pair are replaced together.
        
        Args:
            config: Configuration settings including input/output paths and replacement ratio
            
        Raises:
            GeometryGenerationError: If replacement fails for any reason
        """
        try:
            input_pdb = str(config.replace_file) if config.replace_file else None
            if not input_pdb or not os.path.exists(input_pdb):
                raise GeometryGenerationError(
                    message=f"Input PDB file not found: {input_pdb}",
                    error_code="GEO_ERR_004"
                )
                
            output_pdb = f"{config.output}.pdb"
            system_type = "NC"
            
            LOG.info(f"- Input PDB: {input_pdb}")

            model_count = self._split_pdb_into_models(input_pdb, system_type)
            
            temp_crosslinks_file = os.path.join(system_type, "temp_crosslinks.txt")
            with open(temp_crosslinks_file, 'w') as f:
                for model_id in range(model_count):
                    model_path = os.path.join(system_type, f"{model_id}.caps.pdb")
                    if not os.path.exists(model_path):
                        continue
                    
                    LOG.debug(f"Processing model file: {model_path}")
                    atom_count = 0
                    crosslink_count = 0
                    
                    with open(model_path, 'r') as model_file:
                        for line in model_file:
                            if line.startswith(("ATOM", "HETATM")):
                                atom_count += 1
                                
                                if len(line) >= 20:
                                    resname = line[17:20].strip()
                                    
                                    if resname in ['L4Y', 'L5Y', 'L4X', 'L5X', 'LY4', 'LY5', 'LX4', 'LX5',
                                                'LGX', 'LPS', 'AGS', 'APD', 'DHL', 'HYL']:
                                        f.write(f"{model_id}|{line}")
                                        crosslink_count += 1
                                        
                    LOG.debug(f"Model {model_id}: {atom_count} atoms total, {crosslink_count} crosslink atoms found")
            
            crosslink_residues = {}
            with open(temp_crosslinks_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('|', 1)
                    if len(parts) != 2:
                        continue
                        
                    model_id = int(parts[0])
                    atom_line = parts[1]
                    
                    atom_name = atom_line[12:16].strip()
                    resname = atom_line[17:20].strip()
                    resid = atom_line[22:26].strip()
                    chain = atom_line[21]
                    x = float(atom_line[30:38])
                    y = float(atom_line[38:46])
                    z = float(atom_line[46:54])
                    
                    key = (model_id, resname, resid, chain)
                    
                    if key not in crosslink_residues:
                        crosslink_residues[key] = {
                            'model_id': model_id,
                            'resname': resname,
                            'resid': resid,
                            'chain': chain,
                            'atoms': [],
                            'position': None
                        }
                    
                    crosslink_residues[key]['atoms'].append({
                        'atom_name': atom_name,
                        'position': [x, y, z]
                    })
            
            for key, residue in crosslink_residues.items():
                atoms = residue['atoms']
                if not atoms:
                    continue
                    
                x_sum = sum(atom['position'][0] for atom in atoms)
                y_sum = sum(atom['position'][1] for atom in atoms)
                z_sum = sum(atom['position'][2] for atom in atoms)
                
                residue['position'] = [
                    x_sum / len(atoms),
                    y_sum / len(atoms),
                    z_sum / len(atoms)
                ]
            
            crosslink_residues = {k: v for k, v in crosslink_residues.items() if v['position'] is not None}
            
            pair_types = [
                (['L4Y', 'L4X', 'LY4', 'LX4'], ['L5Y', 'L5X', 'LY5', 'LX5']),
                (['LGX', 'LPS'], ['AGS', 'APD']),
                (['DHL'], ['HYL'])
            ]
            
            residues_by_type = {}
            for key, residue in crosslink_residues.items():
                resname = residue['resname']
                if resname not in residues_by_type:
                    residues_by_type[resname] = []
                residues_by_type[resname].append((key, residue))
            
            LOG.debug(f"Found residue types: {list(residues_by_type.keys())}")
            
            all_pairs = []
            
            for type1_list, type2_list in pair_types:
                type1_residues = []
                for resname in type1_list:
                    if resname in residues_by_type:
                        type1_residues.extend(residues_by_type[resname])
                
                type2_residues = []
                for resname in type2_list:
                    if resname in residues_by_type:
                        type2_residues.extend(residues_by_type[resname])
                
                LOG.debug(f"Looking for pairs between {type1_list} ({len(type1_residues)}) and {type2_list} ({len(type2_residues)})")
                
                if not type1_residues or not type2_residues:
                    continue

                for i, (key1, residue1) in enumerate(type1_residues):
                    for j, (key2, residue2) in enumerate(type2_residues):
                        try:
                            pos1 = residue1['position']
                            pos2 = residue2['position']
                            dist = self._calculate_distance(pos1, pos2)
                            
                            if dist <= 10.0:
                                pair = {
                                    'members': [residue1, residue2],
                                    'distance': dist,
                                    'keys': [key1, key2]
                                }
                                all_pairs.append(pair)
                                LOG.debug(f"Found pair at distance {dist:.2f}Å: "
                                        f"{residue1['resname']} {residue1['resid']}{residue1['chain']} (model {residue1['model_id']}) - "
                                        f"{residue2['resname']} {residue2['resid']}{residue2['chain']} (model {residue2['model_id']})")
                        except Exception as e:
                            LOG.warning(f"Error calculating distance between {residue1['resname']} {residue1['resid']}{residue1['chain']} and "
                                    f"{residue2['resname']} {residue2['resid']}{residue2['chain']}: {e}")
            
            all_pairs.sort(key=lambda p: p['distance'])
            
            used_residues = set()
            unique_pairs = []
            
            for pair in all_pairs:
                key1 = pair['keys'][0]
                key2 = pair['keys'][1]
                
                if key1 not in used_residues and key2 not in used_residues:
                    unique_pairs.append(pair)
                    used_residues.add(key1)
                    used_residues.add(key2)
            
            num_to_replace = min(
                math.ceil(len(unique_pairs) * config.ratio_replace / 100),
                len(unique_pairs)
            )
            
            LOG.debug(f"Will replace {num_to_replace} pairs out of {len(unique_pairs)} (ratio: {config.ratio_replace}%)")
            
            pairs_to_replace = random.sample(unique_pairs, num_to_replace) if num_to_replace > 0 else []
            
            instructions = []
            
            for pair in pairs_to_replace:
                for residue in pair['members']:
                    model_id = residue['model_id']
                    instructions.append(f"{model_id}.caps.pdb LYS {residue['resid']} {residue['chain']}")
                    LOG.debug(f"Will replace {residue['resname']} {residue['resid']}{residue['chain']} in model {model_id}")
            
            if not instructions:
                LOG.warning("No replacement instructions generated")
                LOG.info("Creating output file without modifications since no replacements were made")
                with open(output_pdb, 'w') as out:
                    with open(input_pdb, 'r') as in_file:
                        first_line = in_file.readline().strip()
                        if first_line.startswith("CRYST1"):
                            out.write(first_line + "\n")
                        else:
                            in_file.seek(0)
                            out.write("REMARK   Generated by colbuilder - no replacements made\n")
                        
                        for line in in_file:
                            if line.startswith(("ATOM", "HETATM", "TER")):
                                out.write(line)
                    
                    out.write("END\n")

                return
            
            success = self._run_chimera_replacement(
                instructions=instructions,
                chimera_scripts_dir=config.CHIMERA_SCRIPTS_DIR,
                system_type=system_type,
                output_pdb=None
            )
            
            if success:
                try:
                    from .system import System
                    from .crystal import Crystal
                    
                    crystal = Crystal(pdb=input_pdb)
                    system = System(crystal=crystal)
                    
                    pdb_files = [f for f in os.listdir(system_type) if f.endswith('.pdb') and f != "temp_crosslinks.pdb"]
                    pdb_files.sort(key=lambda x: int(x.split('.')[0]))
                    
                    with open(output_pdb, 'w') as out:
                        with open(input_pdb, 'r') as in_file:
                            first_line = in_file.readline().strip()
                            if first_line.startswith("CRYST1"):
                                out.write(first_line + "\n")
                            else:
                                out.write("REMARK   Generated by colbuilder direct replacement\n")
                        
                        for pdb_file in pdb_files:
                            pdb_path = os.path.join(system_type, pdb_file)
                            with open(pdb_path, 'r') as f:
                                for line in f:
                                    if line.startswith(("ATOM", "HETATM", "TER")):
                                        if line.startswith('HETATM'):
                                            line = 'ATOM  ' + line[6:]
                                        out.write(line)
                        
                        out.write("END\n")
                    
                except Exception as e:
                    LOG.error(f"Error creating output PDB with System class: {str(e)}")
                    LOG.info("Falling back to direct file combining method")
                    
                    with open(output_pdb, 'w') as out:
                        out.write("REMARK   Generated by colbuilder direct replacement\n")
                        try:
                            with open(input_pdb, 'r') as in_file:
                                for line in in_file:
                                    if line.startswith("CRYST1"):
                                        out.write(line)
                                        break
                        except Exception as e:
                            LOG.warning(f"Could not retrieve CRYST1 info: {e}")
                        
                        pdb_files = [f for f in os.listdir(system_type) if f.endswith('.pdb') and f != "temp_crosslinks.pdb"]
                        pdb_files.sort(key=lambda x: int(x.split('.')[0]))
                        
                        for pdb_file in pdb_files:
                            pdb_path = os.path.join(system_type, pdb_file)
                            with open(pdb_path, 'r') as f:
                                for line in f:
                                    if line.startswith(("ATOM", "HETATM", "TER")):
                                        if line.startswith('HETATM'):
                                            line = 'ATOM  ' + line[6:]
                                        out.write(line)
                        
                        out.write("END\n")
                        
                    LOG.info(f"Created output PDB file using fallback method: {output_pdb}")
            else:
                LOG.error("Chimera replacement failed")
                raise GeometryGenerationError(
                    message="Chimera replacement process failed",
                    error_code="GEO_ERR_004"
                )
        
        except GeometryGenerationError:
            raise
        except Exception as e:
            LOG.error(f"Unexpected error in direct replacement: {str(e)}")
            traceback.print_exc()
            raise GeometryGenerationError(
                message=f"Unexpected error in direct replacement: {str(e)}",
                error_code="GEO_ERR_004"
            )
                
    def _calculate_fibril_bounds(self, system: Any, fibril_length: float) -> Tuple[float, float]:
        """
        Calculate the z-coordinate bounds of the collagen fibril.
        
        Determines an appropriate z-coordinate range to focus crosslink replacement,
        based on model centers of geometry and the specified fibril length.
        
        Args:
            system: System containing models with spatial information
            fibril_length: Length of the fibril in nanometers
            
        Returns:
            Tuple of (z_min, z_max) as floating point values defining the bounds
        """
        try:
            z_values = []
            
            for model_id in system.get_models():
                model = system.get_model(model_id=model_id)

                try:
                    cog = model.get_cog()
                    z_center = self._get_z_position(cog)
                    z_values.append(z_center)
                except Exception as e:
                    LOG.warning(f"Could not get COG for model {model_id}: {e}")
            
            if not z_values:
                LOG.warning("No valid z-coordinates found in system models or crosslinks")
                return (0.0, 10000.0)
                
            min_z = min(z_values) - 500.0
            max_z = max(z_values) + 500.0
            
            if max_z - min_z > 5000.0:
                model = system.get_model(model_id=0.0)
                cog = model.get_cog()
                z_center = self._get_z_position(cog)
                
                fibril_length_angstroms = fibril_length * 10.0
                
                z_min = z_center - fibril_length_angstroms * 5
                z_max = z_center + fibril_length_angstroms * 5
            else:
                z_min = min_z
                z_max = max_z
                
            return (z_min, z_max)
        except Exception as e:
            LOG.error(f"Error calculating fibril bounds: {e}")
            return (0.0, 10000.0)
          
    def _get_z_position(self, position: Any) -> float:
        """
        Extract the z-component from a position vector.
        
        Safely extracts the z-coordinate from various position representations
        (numpy arrays, lists, tuples) and converts it to a Python float.
        
        Args:
            position: Position object which may be a numpy array, list, tuple or similar
            
        Returns:
            Z-coordinate as a Python float
            
        Raises:
            IndexError, TypeError, ValueError: If extraction fails
        """
        try:
            if hasattr(position, 'shape'):
                if len(position.shape) == 1 and position.shape[0] >= 3:
                    return float(position[2])
                elif len(position.shape) == 2 and position.shape[1] >= 3:
                    return float(position[0][2])
            
            if hasattr(position, '__getitem__'):
                z_val = position[2]
                if hasattr(z_val, '__getitem__') and hasattr(z_val, '__len__') and len(z_val) > 0:
                    return float(z_val[0])
                return float(z_val)
            
            return float(position)
            
        except (IndexError, TypeError, ValueError) as e:
            LOG.error(f"Error extracting z position: {e}, position type: {type(position)}, value: {position}")
            raise
    
    def _calculate_distance(self, pos1: Any, pos2: Any) -> float:
        """
        Calculate Euclidean distance between two positions.
        
        Handles various position representations (numpy arrays, lists, tuples)
        and computes the 3D Euclidean distance between them.
        
        Args:
            pos1: First position as a vector-like object
            pos2: Second position as a vector-like object
            
        Returns:
            Distance as a Python float
            
        Raises:
            IndexError, TypeError, ValueError: If position extraction fails
        """
        try:
            x1 = self._get_component(pos1, 0)
            y1 = self._get_component(pos1, 1)
            z1 = self._get_component(pos1, 2)
            
            x2 = self._get_component(pos2, 0)
            y2 = self._get_component(pos2, 1)
            z2 = self._get_component(pos2, 2)
            
            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1
            
            return (dx*dx + dy*dy + dz*dz)**0.5
            
        except (IndexError, TypeError, ValueError) as e:
            LOG.error(f"Error calculating distance: {e}")
            LOG.error(f"pos1 type: {type(pos1)}, value: {pos1}")
            LOG.error(f"pos2 type: {type(pos2)}, value: {pos2}")
            raise
    
    def _get_component(self, position: Any, index: int) -> float:
        """
        Extract a specific component from a position vector.
        
        Safely extracts the specified component (x=0, y=1, z=2) from various 
        position representations and converts it to a Python float.
        
        Args:
            position: Position object which may be a numpy array, list, tuple or similar
            index: Component index to extract (0=x, 1=y, 2=z)
            
        Returns:
            Component value as a Python float
            
        Raises:
            IndexError, TypeError, ValueError: If extraction fails
        """
        try:
            if hasattr(position, 'shape'):
                if len(position.shape) == 1 and position.shape[0] > index:
                    return float(position[index])
                elif len(position.shape) == 2 and position.shape[1] > index:
                    return float(position[0][index])
            
            if hasattr(position, '__getitem__'):
                val = position[index]
                if hasattr(val, '__getitem__') and hasattr(val, '__len__') and len(val) > 0:
                    return float(val[0])
                return float(val)
            
            return float(position)
            
        except (IndexError, TypeError, ValueError) as e:
            LOG.error(f"Error extracting component {index}: {e}")
            LOG.error(f"position type: {type(position)}, value: {position}")
            raise
    
    def _write_replace_file(self, system: Any, replace_file: str) -> None:
        """
        Write replacement instructions to a file for Chimera.
        
        Creates a file containing replacement instructions for crosslinks
        marked with state='replace' in the system. Each line specifies a model,
        residue name, residue ID, and chain for Chimera to process.
        
        Args:
            system: System with marked crosslinks
            replace_file: Path to the output file for instructions
            
        Raises:
            GeometryGenerationError: If the file cannot be written
        """
        try:
            with open(replace_file, 'w') as f:
                for key in system.get_models():
                    model = system.get_model(model_id=key)
                    if not hasattr(model, 'crosslink') or not model.crosslink:
                        continue
                        
                    for cross in model.crosslink:
                        if cross.state == 'replace':
                            f.write(f"{int(key)}.caps.pdb {cross.resname} {cross.resid} {cross.chain}\n")
        except Exception as e:
            LOG.error(f"Failed to write replacement file: {e}")
            raise GeometryGenerationError(
                message=f"Failed to write replacement file: {str(e)}",
                error_code="GEO_ERR_004"
            )
    
    async def _run_replacement_with_chimera(
        self, 
        system: Any, 
        config: ColbuilderConfig, 
        replace_file: str
    ) -> None:
        """
        Execute Chimera to perform crosslink replacements.
        
        Runs Chimera with appropriate scripts to replace crosslinks with lysines
        according to the instructions in the replacement file.
        
        Args:
            system: System with crosslinks to replace
            config: Configuration settings for Chimera execution
            replace_file: Path to file with replacement instructions
            
        Raises:
            GeometryGenerationError: If Chimera execution fails
        """
        try:
            LOG.debug("Running Chimera to replace crosslinks with lysines")
            
            model_zero = system.get_model(model_id=0.0)
            if not hasattr(model_zero, 'type'):
                raise GeometryGenerationError(
                    message="Model does not have a 'type' attribute needed for replacement",
                    error_code="GEO_ERR_004"
                )
                
            system_type = model_zero.type
            
            os.makedirs(system_type, exist_ok=True)
            
            base_file = os.path.splitext(replace_file)[0]
            success = self._run_chimera_replacement(
                instructions=self._read_replacement_instructions(replace_file),
                chimera_scripts_dir=config.CHIMERA_SCRIPTS_DIR,
                system_type=system_type
            )
            
            if not success:
                raise GeometryGenerationError(
                    message="Chimera replacement failed",
                    error_code="GEO_ERR_004"
                )
                
            LOG.debug("Chimera replacement completed successfully")
            
        except subprocess.SubprocessError as e:
            LOG.error(f"Subprocess error during Chimera execution: {str(e)}")
            raise GeometryGenerationError(
                message=f"Chimera execution failed: {str(e)}",
                error_code="GEO_ERR_004"
            )
    
    def _read_replacement_instructions(self, replace_file: str) -> List[str]:
        """
        Read replacement instructions from a file.
        
        Parses the replacement instructions file and returns a list of instruction strings.
        
        Args:
            replace_file: Path to the replacement instructions file
            
        Returns:
            List of replacement instruction strings
            
        Raises:
            GeometryGenerationError: If the file cannot be read
        """
        instructions = []
        try:
            with open(replace_file, 'r') as f:
                instructions = [line.strip() for line in f if line.strip()]
        except Exception as e:
            LOG.error(f"Error reading replacement file: {e}")
            raise GeometryGenerationError(
                message=f"Failed to read replacement file: {str(e)}",
                error_code="GEO_ERR_004"
            )
        return instructions
    
    def _run_chimera_replacement(
        self,
        instructions: List[str],
        chimera_scripts_dir: str,
        system_type: str,
        output_pdb: Optional[str] = None
    ) -> bool:
        """
        Run Chimera to perform replacements on PDB files.
        
        Executes Chimera with the swapaa.py script to replace residues in the
        specified PDB files according to the provided instructions.
        
        Args:
            instructions: List of replacement instructions
            chimera_scripts_dir: Directory containing Chimera scripts
            system_type: Type of system (directory name)
            output_pdb: Optional path for output PDB file
            
        Returns:
            True if successful, False if replacement failed
        """
        try:
            replace_file = "replace.txt"
            with open(replace_file, 'w') as f:
                f.write("\n".join(instructions))
                
            LOG.debug(f"Created replacement file: {replace_file}")
            
            LOG.debug("Running Chimera to replace crosslinks")
            swapaa_script = os.path.join(chimera_scripts_dir, "swapaa.py")
            
            if not os.path.exists(swapaa_script):
                raise GeometryGenerationError(
                    message=f"Chimera swapaa script not found: {swapaa_script}",
                    error_code="GEO_ERR_004"
                )
                
            base_file = replace_file[:-4] if replace_file.endswith('.txt') else replace_file
            
            cmd = f"chimera --nogui --silent --script \"{swapaa_script} {base_file} {system_type}\""
            
            LOG.debug(f"Running command: {cmd}")
            
            result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            LOG.debug(f"Chimera returned: {result.returncode}")
            
            if result.returncode != 0:
                LOG.error(f"Chimera failed with return code {result.returncode}")
                if result.stderr:
                    LOG.error(f"Chimera stderr: {result.stderr}")
                if result.stdout:
                    LOG.info(f"Chimera stdout: {result.stdout}")
                    
                return False
                
            if output_pdb:
                with open(output_pdb, 'w') as out:
                    out.write("REMARK   Generated by colbuilder replacement\n")
                    
                    pdb_files = [f for f in os.listdir(system_type) if f.endswith('.pdb')]
                    pdb_files.sort(key=lambda x: int(x.split('.')[0]))
                    
                    for pdb_file in pdb_files:
                        pdb_path = os.path.join(system_type, pdb_file)
                        with open(pdb_path, 'r') as f:
                            for line in f:
                                if line.startswith(("ATOM", "HETATM", "TER")):
                                    out.write(line)
                    
                    out.write("END\n")
                
            LOG.debug("Chimera replacement completed successfully")
            return True
            
        except Exception as e:
            LOG.error(f"Error during Chimera execution: {str(e)}")
            traceback.print_exc()
            return False
    
    def _split_pdb_into_models(self, input_pdb: str, system_type: str) -> int:
        """
        Split a PDB file into individual triple helix models.
        
        Divides the PDB file into individual models using the repeating A-B-C
        chain pattern that is typical of collagen triple helices.
        
        Args:
            input_pdb: Path to input PDB file
            system_type: Type of system/directory to create
            
        Returns:
            Number of models created
            
        Raises:
            GeometryGenerationError: If the input PDB file is not a valid Colbuilder structure
        """
        is_colbuilder_pdb = False
        with open(input_pdb, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith("CRYST1"):
                is_colbuilder_pdb = True
                
        if not is_colbuilder_pdb:
            raise GeometryGenerationError(
                message="Input PDB file does not appear to be a valid Colbuilder-generated structure. " 
                        "Direct replacement only works on PDB files generated by the Colbuilder geometry module.",
                error_code="GEO_ERR_004"
            )
        
        with open(input_pdb, 'r') as f:
            all_lines = f.readlines()
        
        cryst_line = None
        for line in all_lines:
            if line.startswith("CRYST1"):
                cryst_line = line
                break
        
        atom_lines = [line for line in all_lines if line.startswith(("ATOM", "HETATM", "TER"))]
        
        if not atom_lines:
            raise GeometryGenerationError(
                message="Input PDB file does not contain any ATOM, HETATM, or TER records.",
                error_code="GEO_ERR_004"
            )
        
        models = []
        current_model = []
        chain_sequence = []
        
        for line in atom_lines:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                
                if not chain_sequence or chain != chain_sequence[-1]:
                    chain_sequence.append(chain)
                    
                    if len(chain_sequence) > 3 and chain_sequence[-4:] == ['A', 'B', 'C', 'A']:
                        models.append(current_model)
                        current_model = []
                        chain_sequence = ['A']
            
            current_model.append(line)
        
        if current_model:
            models.append(current_model)
        
        if not models:
            LOG.warning("Could not split into multiple models. Using entire PDB as one model.")
            models = [atom_lines]
        
        LOG.debug(f"Split PDB into {len(models)} models")
        
        if not os.path.exists(system_type):
            os.makedirs(system_type)
        
        for i, model_lines in enumerate(models):
            caps_file = os.path.join(system_type, f"{i}.caps.pdb")
            with open(caps_file, 'w') as f:
                if cryst_line:
                    f.write(cryst_line)
                
                f.writelines(model_lines)
                
                if not model_lines[-1].startswith("TER"):
                    f.write("TER\n")
        
        return len(models)

    def _select_replacements(self, system: Any, ratio_replace: float, z_bounds: Tuple[float, float]) -> None:
        """
        Select crosslinks to replace based on the given ratio and spatial constraints.
        
        Identifies complete crosslink pairs within the specified z-bounds, ensures
        they are properly paired, and marks a subset for replacement based on the
        specified ratio.
        
        Args:
            system: System containing models and crosslinks
            ratio_replace: Percentage of crosslinks to replace (0-100)
            z_bounds: Tuple of (z_min, z_max) defining the region to focus on
        """
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            if hasattr(model, 'crosslink') and model.crosslink:
                for crosslink in model.crosslink:
                    crosslink.state = 'none'
        
        total_system_crosslinks = 0
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            if hasattr(model, 'crosslink') and model.crosslink:
                total_system_crosslinks += len(model.crosslink)
        
        LOG.debug(f"Total crosslinks in system: {total_system_crosslinks}")
        
        crosslinks_in_bounds = []
        for model_id in system.get_models():
            model = system.get_model(model_id=model_id)
            if not hasattr(model, 'crosslink') or not model.crosslink:
                continue
                
            for crosslink in model.crosslink:
                try:
                    z_pos = self._get_z_position(crosslink.position)
                    if z_bounds[0] <= z_pos <= z_bounds[1]:
                        crosslinks_in_bounds.append({
                            'model_id': model_id,
                            'model': model,
                            'crosslink': crosslink,
                            'type': getattr(crosslink, 'type', 'D')
                        })
                except Exception as e:
                    LOG.warning(f"Error checking crosslink bounds: {e}")
        
        crosslinks_by_type = {}
        for cross_ref in crosslinks_in_bounds:
            crosslink = cross_ref['crosslink']
            resname = crosslink.resname
            if resname not in crosslinks_by_type:
                crosslinks_by_type[resname] = []
            crosslinks_by_type[resname].append(cross_ref)
        
        LOG.debug(f"Found crosslink types: {list(crosslinks_by_type.keys())}")
        
        pair_types = [
            (['L4Y', 'L4X', 'LY4', 'LX4'], ['L5Y', 'L5X', 'LY5', 'LX5']),
            (['LGX', 'LPS'], ['AGS', 'APD']),
            (['DHL'], ['HYL'])
        ]
        
        pairs = []
        
        for type1_list, type2_list in pair_types:
            type1_crosslinks = []
            for resname in type1_list:
                if resname in crosslinks_by_type:
                    type1_crosslinks.extend(crosslinks_by_type[resname])
            
            type2_crosslinks = []
            for resname in type2_list:
                if resname in crosslinks_by_type:
                    type2_crosslinks.extend(crosslinks_by_type[resname])
            
            if not type1_crosslinks or not type2_crosslinks:
                continue
            
            used_type1 = set()
            used_type2 = set()
            
            for i, cross_ref1 in enumerate(type1_crosslinks):
                if i in used_type1:
                    continue
                    
                crosslink1 = cross_ref1['crosslink']
                min_dist = float('inf')
                closest_idx = -1
                
                for j, cross_ref2 in enumerate(type2_crosslinks):
                    if j in used_type2:
                        continue
                        
                    crosslink2 = cross_ref2['crosslink']
                    
                    try:
                        dist = self._calculate_distance(crosslink1.position, crosslink2.position)
                        if dist < min_dist and dist <= 5.0:
                            min_dist = dist
                            closest_idx = j
                    except Exception as e:
                        LOG.warning(f"Error calculating distance: {e}")
                
                if closest_idx != -1:
                    pairs.append([cross_ref1, type2_crosslinks[closest_idx]])
                    used_type1.add(i)
                    used_type2.add(closest_idx)
                    
                    LOG.debug(f"Found pair: {crosslink1.resname} {crosslink1.resid}{crosslink1.chain} + "
                            f"{type2_crosslinks[closest_idx]['crosslink'].resname} {type2_crosslinks[closest_idx]['crosslink'].resid}{type2_crosslinks[closest_idx]['crosslink'].chain} "
                            f"(distance: {min_dist:.2f}Å)")
        
        num_to_replace = min(
            math.ceil(len(pairs) * ratio_replace / 100),
            len(pairs)
        )
        
        LOG.debug(f"Will replace {num_to_replace} pairs out of {len(pairs)} (ratio: {ratio_replace}%)")
        
        random.shuffle(pairs)
        pairs_to_replace = pairs[:num_to_replace]
        
        replaced_count = 0
        for pair in pairs_to_replace:
            for cross_ref in pair:
                cross_ref['crosslink'].state = 'replace'
                replaced_count += 1
                
                model_id = cross_ref['model_id']
                crosslink = cross_ref['crosslink']
                LOG.debug(f"Marked {crosslink.resname} {crosslink.resid}{crosslink.chain} in model {model_id} for replacement")
        
        for pair in pairs_to_replace:
            for cross_ref in pair:
                crosslink = cross_ref['crosslink']
                model = cross_ref['model']
                
                for other in model.crosslink:
                    if other.state != 'none' or other in [p['crosslink'] for p in pair]:
                        continue
                        
                    try:
                        distance = self._calculate_distance(crosslink.position, other.position)
                        if 5.0 < distance <= 15.0:
                            other.state = 'protect'
                    except Exception as e:
                        LOG.warning(f"Error protecting nearby crosslink: {e}")
        
        if total_system_crosslinks > 0:
            achieved_ratio = (replaced_count / total_system_crosslinks) * 100
            LOG.debug(f"Achieved replacement ratio: {achieved_ratio:.2f}% (target: {ratio_replace}%)")

# Backward compatibility functions

async def replace_in_system(system: Any, config: ColbuilderConfig) -> Any:
    """
    Backward compatibility function for system-based crosslink replacement.
    
    This is a standalone function that creates a CrosslinkReplacer and calls
    its replace_in_system method, maintaining backward compatibility with
    code that uses the function-based API.
    
    Args:
        system: System containing the models to be modified
        config: Configuration for replacement
            
    Returns:
        Modified system with replaced crosslinks
    """
    replacer = CrosslinkReplacer()
    return await replacer.replace_in_system(system, config)

async def direct_replace_geometry(config: ColbuilderConfig) -> None:
    """
    Backward compatibility function for direct PDB-based crosslink replacement.
    
    This is a standalone function that creates a CrosslinkReplacer and calls
    its replace_direct method, maintaining backward compatibility with
    code that uses the function-based API.
    
    Args:
        config: Configuration settings for the replacement
    """
    replacer = CrosslinkReplacer()
    await replacer.replace_direct(config)