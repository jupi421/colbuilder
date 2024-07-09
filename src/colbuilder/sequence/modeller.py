import os
import sys
from modeller import *
from modeller.automodel import *
from modeller import Selection

class Modeller:
    def __init__(self, fasta_file, aligned_sequences, file_prefix, original_dir):
        self.fasta_file = fasta_file
        self.aligned_sequences = aligned_sequences
        self.file_prefix = file_prefix
        self.original_dir = original_dir
        self.modeller_pdb = ''
        modeller.log.minimal()

    def read_sequences_from_fasta(self, filename):
        """Read sequences from a FASTA file and return as a dictionary."""
        sequences = {}
        with open(filename) as f:
            for record in SeqIO.parse(f, "fasta"):
                sequences[record.id] = str(record.seq)
        if not sequences:
            raise ValueError(f"FASTA file '{filename}' is empty or incorrectly formatted.")
        return sequences

    def get_template_pdb_file(self):
        """Fetch the template.pdb file."""
        template_pdb = os.path.join(os.path.dirname(self.file_prefix), "template.pdb")
        if not os.path.exists(template_pdb):
            raise FileNotFoundError(f"template.pdb not found in the directory: {os.path.dirname(self.file_prefix)}")
        return template_pdb

    def extract_last_atom_serial_number(self, pdb_file):
        """Extracts the serial number from the last ATOM line in a PDB file, before the final TER line."""
        last_atom_serial_number = None
        with open(pdb_file, "r") as file:
            for line in reversed(file.readlines()):
                if line.startswith("ATOM"):
                    last_atom_serial_number = line[6:11].strip()
                    break
        return last_atom_serial_number
        
    def write_modeller_formatted_output(self, fasta_sequences, msa_sequences, output_file, template_pdb, end_atom):
        """Write output in the format required by Modeller."""
        pdb_file_name_without_ext = os.path.splitext(os.path.basename(template_pdb))[0]
        
        with open(output_file, "w") as out_file:
            # Writing the template sequence
            out_file.write(">P1;template\n")
            out_file.write(f"structure:{pdb_file_name_without_ext}:1:A:+{end_atom}:C:::-1.00:-1.00\n")
            out_file.write("/".join(fasta_sequences.values()) + "*\n")

            # Writing the target sequence
            out_file.write("\n>P1;target\n")
            out_file.write(f"sequence:{self.collagen_type}-{self.register}:: :: :::-1.00:-1.00\n")
            out_file.write("/".join(msa_sequences.values()) + "*\n")
            
    def prepare_alignment(self):
        """
        Prepare the alignment from aligned sequences as input for modeller
        """
        fasta_sequences = self.read_sequences_from_fasta(self.fasta_file)
        msa_sequences = self.read_sequences_from_fasta(self.aligned_sequences)
        
        template_pdb = self.get_template_pdb_file()
        end_atom = self.extract_last_atom_serial_number(template_pdb)
        
        output_file = f"{self.file_prefix}_modeller.ali"
        self.write_modeller_formatted_output(fasta_sequences, msa_sequences, output_file, template_pdb, end_atom)

    def run_modeller(self, alnfile, knowns, sequence, model_dir, out_model_file):
        log.verbose()
        env = Environ(
            rand_seed=-8123,
            restyp_lib_file=f"{self.original_dir}/parameters/restyp_mod.lib",
            copy=None,
        )
        
        env.io.atom_files_directory = ["."]
        env.io.hetatm = True

        env.libs.topology.read(f"{self.original_dir}/parameters/top_heav_mod.lib")
        env.libs.parameters.read(f"{self.original_dir}/parameters/par_mod.lib")

        a = MyModel(
            env,
            alnfile=alnfile,
            knowns=knowns,
            sequence=sequence
        )

        a.very_fast()
        a.starting_model = 1
        a.ending_model = 1

        os.chdir(model_dir)

        a.make()
        a.write(file=out_model_file)

        self.modeller_pdb = out_model_file

    def execute_modeller(self):
        alnfile = f"{self.file_prefix}_modeller.ali"
        knowns = "template"
        sequence = "target"
        model_dir = os.path.dirname(self.file_prefix)
        out_model_file = f"{self.file_prefix}_final_model.pdb"

        self.run_modeller(
            alnfile,
            knowns,
            sequence,
            model_dir,
            out_model_file
        )
