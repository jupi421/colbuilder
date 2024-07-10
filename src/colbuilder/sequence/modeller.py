import os
from modeller import *
from modeller.automodel import *
from Bio import SeqIO

class Modeller:
    def __init__(self, aligned_fasta, output_prefix, template_pdb, restyp_lib, top_heav_lib, par_mod_lib):
        """
        Initialize the Modeller class.
        
        :param aligned_fasta: Path to the FASTA file containing aligned sequences
        :param output_prefix: Prefix for output files (without extension)
        :param template_pdb: Path to the template PDB file
        :param restyp_lib: Path to the residue type library file for MODELLER
        :param top_heav_lib: Path to the topology library file for MODELLER
        :param par_mod_lib: Path to the parameter library file for MODELLER
        """
        self.aligned_fasta = aligned_fasta
        self.output_prefix = output_prefix
        self.template_pdb = template_pdb
        self.restyp_lib = restyp_lib
        self.top_heav_lib = top_heav_lib
        self.par_mod_lib = par_mod_lib
        self.output_pdb = None

    def read_sequences_from_fasta(self, filename):
        """Read sequences from a FASTA file and return as a dictionary."""
        sequences = {}
        with open(filename) as f:
            for record in SeqIO.parse(f, "fasta"):
                sequences[record.id] = str(record.seq)
        if not sequences:
            raise ValueError(f"FASTA file '{filename}' is empty or incorrectly formatted.")
        return sequences

    def extract_last_atom_serial_number(self, pdb_file):
        """Extracts the serial number from the last ATOM line in a PDB file."""
        with open(pdb_file, "r") as file:
            for line in reversed(file.readlines()):
                if line.startswith("ATOM"):
                    return line[6:11].strip()
        return None

    def write_modeller_formatted_output(self, sequences, output_file, end_atom):
        """Write output in the format required by Modeller."""
        pdb_file_name_without_ext = os.path.splitext(os.path.basename(self.template_pdb))[0]
       
        with open(output_file, "w") as out_file:
            # Writing the template sequence
            out_file.write(">P1;template\n")
            out_file.write(f"structure:{pdb_file_name_without_ext}:1:A:+{end_atom}:C:::-1.00:-1.00\n")
            out_file.write("/".join(sequences.values()) + "*\n")

            # Writing the target sequence
            out_file.write("\n>P1;target\n")
            out_file.write(f"sequence:{self.output_prefix}:: :: :::-1.00:-1.00\n")
            out_file.write("/".join(sequences.values()) + "*\n")
           
    def prepare_alignment(self):
        """Prepare the alignment from aligned sequences as input for modeller"""
        sequences = self.read_sequences_from_fasta(self.aligned_fasta)
        end_atom = self.extract_last_atom_serial_number(self.template_pdb)
        output_file = f"{self.output_prefix}_modeller.ali"
        self.write_modeller_formatted_output(sequences, output_file, end_atom)

    def run_modeller(self):
        """Run MODELLER to generate the model"""
        log.verbose()
        env = Environ(
            rand_seed=-8123,
            restyp_lib_file=self.restyp_lib,
            copy=None,
        )
       
        env.io.atom_files_directory = ["."]
        env.io.hetatm = True

        env.libs.topology.read(self.top_heav_lib)
        env.libs.parameters.read(self.par_mod_lib)

        a = AutoModel(
            env,
            alnfile=f"{self.output_prefix}_modeller.ali",
            knowns="template",
            sequence="target"
        )

        a.very_fast()
        a.starting_model = 1
        a.ending_model = 1

        a.make()
        
        self.output_pdb = f"{self.output_prefix}_final_model.pdb"
        os.rename(a.outputs[0]['name'], self.output_pdb)

    def execute_modeller(self):
        """Execute the entire MODELLER process"""
        self.prepare_alignment()
        self.run_modeller()
