# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0
 
import os
import subprocess
import logging
import shutil

LOG = logging.getLogger(__name__)

class Amber:
    """
    A class to generate topology for the AMBER99 force field.

    This class provides functionality to merge PDB files, write ITP (Include Topology) files,
    create topology files, and generate GRO (Gromos87) files for use with the AMBER99 force field.

    Attributes:
        system: The molecular system being processed.
        ff (str): The force field name, with '.ff' appended.
        is_line (tuple): Tuple of strings representing valid line starts in PDB files.
    """

    def __init__(self, system=None, ff=None):
        """
        Initialize the Amber object.

        Args:
            system: The molecular system to process.
            ff (str): The force field name (without '.ff').
        """
        self.system = system
        self.ff = ff + '.ff'
        self.is_line = ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
    
    def merge_pdbs(self, connect_id=None):
        """
        Merge PDB files according to the connect_id in the system.

        This method combines multiple PDB files associated with a given connect_id
        into a single merged PDB file.

        Args:
            connect_id: The identifier for the connection in the system.

        Returns:
            str or None: The type of the model if successful, None otherwise.

        Raises:
            FileNotFoundError: If an input PDB file is not found.
        """
        LOG.debug(f"Merging PDFs for connect_id: {connect_id}")
        model = self.system.get_model(model_id=connect_id)
        if model is None or model.connect is None:
            LOG.warning(f"No model or connections found for connect_id: {connect_id}")
            return None

        type_ = model.type
        if not type_:
            LOG.error(f"Type is None for model with connect_id: {connect_id}")
            return None

        os.makedirs(type_, exist_ok=True)
        output_file = os.path.join(type_, f"{int(connect_id)}.merge.pdb")

        if len(model.connect) > 1:
            with open(output_file, 'w') as f:
                for connected_model in model.connect:
                    input_file = os.path.join(type_, f"{int(connected_model)}.caps.pdb")
                    if not os.path.exists(input_file):
                        LOG.error(f"Input file not found: {input_file}")
                        continue
                    with open(input_file, 'r') as infile:
                        f.write("".join(line for line in infile if line.startswith(self.is_line)))
                f.write("END\n")
            LOG.debug(f"Merged PDB written to: {output_file}")
        elif len(model.connect) == 1:  # This is the case for single-model connections
            input_file = os.path.join(type_, f"{int(connect_id)}.caps.pdb")
            if os.path.exists(input_file):
                shutil.copy(input_file, output_file)
            else:
                LOG.error(f"Input file not found for single-model case: {input_file}")
                return None
        else:
            LOG.warning(f"No connections found for connect_id: {connect_id}")
            return None

        return type_
    
    def write_itp(self, itp_file=None):
        """
        Read an ITP file, clean it, and write a new version.

        This method processes an input topology file, removes water topology,
        and writes a cleaned version with only the necessary molecule information.

        Args:
            itp_file (str): Path to the input topology file.

        Raises:
            FileNotFoundError: If the input ITP file is not found.
            PermissionError: If there's no write permission for the output file.
        """
        LOG.debug(f"Writing ITP file: {itp_file}")
        try:
            with open(str(itp_file), 'r') as f:
                itp_model = f.readlines()
        except FileNotFoundError:
            LOG.error(f"Input ITP file not found: {itp_file}")
            raise

        subprocess.run("rm " + str(itp_file), shell=True)
        write = False
        output_file = str(itp_file).replace("top", "itp")
        try:
            with open(output_file, 'w') as f:
                for line in itp_model:
                    if 'Include water topology' in line:
                        break
                    if write:
                        f.write(line)
                    elif 'Protein_chain_A' in line:
                        f.write('[ moleculetype ]\n')
                        f.write(str(itp_file).replace(".top", "") + '  3\n')
                        write = True
            LOG.debug(f"ITP file written: {output_file}")
        except PermissionError:
            LOG.error(f"Permission denied when writing to file: {output_file}")
            raise
    
    def write_topology(self, system=None, topology_file=None, processed_models=None):
        """
        Write a topology file for AMBER99-ILDNP-STAR force field.

        This method generates a comprehensive topology file including
        force field parameters, molecule topologies, and system composition.

        Args:
            system: The molecular system (not used in the current implementation).
            topology_file (str): Path to the output topology file.
            processed_models (list): List of processed model identifiers.

        Raises:
            ValueError: If no processed models are provided.
            PermissionError: If there's no write permission for the output file.
        """
        if not processed_models:
            LOG.error("No processed models to write topology")
            raise ValueError("processed_models cannot be empty")

        LOG.debug(f"Writing topology file: {topology_file}")
        try:
            with open(topology_file, 'w') as f:
                f.write('; Topology for Collagen Microfibril from Colbuilder 2.0\n')
                f.write('#include "./' + self.ff + '/forcefield.itp"\n')
                for model in processed_models:
                    if os.path.exists(f"col_{int(model)}.itp"):
                        f.write(f'#include "col_{int(model)}.itp"\n')
                
                f.write('#include "./' + self.ff + '/ions.itp"\n')
                f.write('#include "./' + self.ff + '/tip3p.itp"\n')
                f.write('\n\n[ system ]\n ;name\nCollagen Microfibril in Water\n\n[ molecules ]\n;name  number\n')
                for model in processed_models:
                    if os.path.exists(f"col_{int(model)}.itp"):
                        f.write(f'col_{int(model)}   1\n')
            LOG.debug(f"Topology file written: {topology_file}")
        except PermissionError:
            LOG.error(f"Permission denied when writing to file: {topology_file}")
            raise

    def write_gro(self, system=None, gro_file=None, processed_models=None):
        """
        Write a GRO (Gromos87) file for the processed models.

        This method combines information from individual GRO files of processed models
        into a single GRO file, updating the total atom count.

        Args:
            system: The molecular system (not used in the current implementation).
            gro_file (str): Path to the output GRO file.
            processed_models (list): List of processed model identifiers.

        Raises:
            ValueError: If no processed models are provided.
            FileNotFoundError: If an input GRO file for a model is not found.
            PermissionError: If there's no write permission for the output file.
        """
        if not processed_models:
            LOG.error("No processed models to write GRO file")
            raise ValueError("processed_models cannot be empty")

        LOG.debug(f"Writing GRO file: {gro_file}")
        gro = []
        total_atoms = 0
        try:
            with open(gro_file, 'w') as f:
                f.write("GROMACS GRO-FILE\n")
                for model in processed_models:
                    model_gro = f"col_{int(model)}.gro"
                    if os.path.exists(model_gro):
                        with open(model_gro, 'r') as model_f:
                            lines = model_f.readlines()
                            total_atoms += int(lines[1])
                            for line in lines[2:-1]:
                                f.write(line)
                        gro = lines  
                        os.remove(model_gro)
                    else:
                        LOG.warning(f"GRO file not found for model: {model}")
                if gro:
                    f.write(gro[-1]) 
            
            with open(gro_file, 'r+') as f:
                content = f.read()
                f.seek(0, 0)
                f.write(f"GROMACS GRO-FILE\n{total_atoms}\n" + content[content.index('\n', content.index('\n') + 1) + 1:])
            
            LOG.debug(f"GRO file written: {gro_file}")
        except PermissionError:
            LOG.error(f"Permission denied when writing to file: {gro_file}")
            raise
        except FileNotFoundError:
            LOG.error(f"One or more input GRO files not found")
            raise