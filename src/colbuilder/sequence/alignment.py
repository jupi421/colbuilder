# alignment.py

import os
import subprocess
import tempfile
import shutil
import logging
from typing import List, Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class Alignment:
    def __init__(self, input_fasta: str, template_fasta: str, output_prefix: str, template_pdb: str):
        self.input_fasta = input_fasta
        self.template_fasta = template_fasta
        self.output_prefix = output_prefix
        self.template_pdb = template_pdb
        self.input_sequences = list(SeqIO.parse(input_fasta, "fasta"))
        self.template_sequences = list(SeqIO.parse(template_fasta, "fasta"))
        self.hydroxyproline_positions: Dict[str, List[int]] = {}
        self.temp_dir = tempfile.mkdtemp()

    def __del__(self):
        shutil.rmtree(self.temp_dir)

    def add_stagger_to_sequences(self, sequences: List[SeqRecord]) -> List[SeqRecord]:
        return [SeqRecord(Seq('-' * (i + 1) + str(seq.seq)), id=seq.id, description=seq.description)
                for i, seq in enumerate(sequences)]

    def process_hydroxyprolines(self, sequence: str) -> Tuple[str, List[int]]:
        seq_list = list(sequence)
        positions = []
        for i, char in enumerate(seq_list):
            if char.upper() == "O":
                positions.append(i)
                seq_list[i] = "P"
        return "".join(seq_list), positions

    def create_position_mapping(self, original_seq: str, aligned_seq: str) -> Dict[int, int]:
        map_original_to_aligned = {}
        aligned_index = 0
        for original_index, char in enumerate(original_seq):
            while aligned_index < len(aligned_seq) and aligned_seq[aligned_index] == "-":
                aligned_index += 1
            if aligned_index < len(aligned_seq):
                map_original_to_aligned[original_index] = aligned_index
                aligned_index += 1
        return map_original_to_aligned

    def restore_hydroxyproline_with_mapping(self, aligned_seq: str, original_positions: List[int], position_mapping: Dict[int, int]) -> str:
        aligned_seq_list = list(aligned_seq)
        for original_pos in original_positions:
            if original_pos in position_mapping:
                new_pos = position_mapping[original_pos]
                aligned_seq_list[new_pos] = 'O'
        return "".join(aligned_seq_list)

    def get_muscle_version(self) -> str:
        try:
            result = subprocess.run(["muscle", "-version"], capture_output=True, text=True, check=True)
            return result.stdout.strip()
        except subprocess.CalledProcessError:
            raise RuntimeError("MUSCLE is not installed or not found in PATH.")

    def align_sequences_with_muscle(self, input_path: str, output_path: str):
        muscle_version = self.get_muscle_version()
        if "3.8" in muscle_version:
            muscle_command = ["muscle", "-in", input_path, "-out", output_path]
        else:
            muscle_command = ["muscle", "-align", input_path, "-output", output_path]
        try:
            subprocess.run(muscle_command, check=True)
            logger.info("Alignment completed.")
        except subprocess.CalledProcessError:
            raise RuntimeError("MUSCLE alignment failed.")

    def extract_last_atom_serial_number(self, pdb_file: str) -> str:
        with open(pdb_file, "r") as file:
            for line in reversed(file.readlines()):
                if line.startswith("ATOM"):
                    return line[6:11].strip()
        return "0"

    def write_modeller_formatted_output(self, sequences: List[SeqRecord], output_file: str, end_atom: str):
        pdb_file_name_without_ext = os.path.splitext(os.path.basename(self.template_pdb))[0]
       
        with open(output_file, "w") as out_file:
            out_file.write(">P1;template\n")
            out_file.write(f"structure:{pdb_file_name_without_ext}:1:A:+{end_atom}:C:::-1.00:-1.00\n")
            out_file.write("|".join(str(seq.seq) for seq in sequences[:3]) + "*\n")

            out_file.write("\n>P1;target\n")
            out_file.write(f"sequence:{self.output_prefix}:: :: :::-1.00:-1.00\n")
            out_file.write("|".join(str(seq.seq) for seq in sequences[3:]) + "*\n")

    def equalize_sequence_lengths(self, sequences: List[SeqRecord]) -> List[SeqRecord]:
        max_length = max(len(seq.seq) for seq in sequences)
        return [SeqRecord(Seq(str(seq.seq).ljust(max_length, '-')), id=seq.id, description=seq.description) for seq in sequences]

    def align_sequences(self) -> Tuple[List[SeqRecord], str, str]:
        # Step 1: Align input sequences
        input_modified_sequences = []
        for seq in self.input_sequences:
            substituted_seq, positions = self.process_hydroxyprolines(str(seq.seq))
            self.hydroxyproline_positions[seq.id] = positions
            modified_seq = SeqRecord(Seq(substituted_seq), id=seq.id, description=seq.description)
            input_modified_sequences.append(modified_seq)

        temp_input_path = os.path.join(self.temp_dir, f"{self.output_prefix}_input.fasta")
        temp_input_aligned_path = os.path.join(self.temp_dir, f"{self.output_prefix}_input_aligned.afa")

        SeqIO.write(input_modified_sequences, temp_input_path, "fasta")
        self.align_sequences_with_muscle(temp_input_path, temp_input_aligned_path)
        
        input_aligned_sequences = list(SeqIO.parse(temp_input_aligned_path, "fasta"))

        # Step 2: Prepare template sequences
        template_modified_sequences = []
        for seq in self.template_sequences:
            substituted_seq, positions = self.process_hydroxyprolines(str(seq.seq))
            self.hydroxyproline_positions[seq.id] = positions
            modified_seq = SeqRecord(Seq(substituted_seq), id=seq.id, description=seq.description)
            template_modified_sequences.append(modified_seq)

        # Step 3: Align input MSA with template
        combined_sequences = template_modified_sequences + input_aligned_sequences
        temp_combined_path = os.path.join(self.temp_dir, f"{self.output_prefix}_combined.fasta")
        temp_final_aligned_path = os.path.join(self.temp_dir, f"{self.output_prefix}_final_aligned.afa")

        SeqIO.write(combined_sequences, temp_combined_path, "fasta")
        self.align_sequences_with_muscle(temp_combined_path, temp_final_aligned_path)

        final_aligned_sequences = list(SeqIO.parse(temp_final_aligned_path, "fasta"))

        # Step 4: Restore hydroxyprolines
        restored_sequences = []
        for seq in final_aligned_sequences:
            if seq.id in self.hydroxyproline_positions:
                original_seq = next(s for s in self.input_sequences + self.template_sequences if s.id == seq.id)
                position_mapping = self.create_position_mapping(str(original_seq.seq), str(seq.seq))
                restored_seq = self.restore_hydroxyproline_with_mapping(str(seq.seq), self.hydroxyproline_positions[seq.id], position_mapping)
                restored_sequences.append(SeqRecord(Seq(restored_seq), id=seq.id, description=seq.description))
            else:
                restored_sequences.append(seq)

        # Step 5: Add stagger to template and input sequences
        template_staggered = self.add_stagger_to_sequences(restored_sequences[:3])
        input_staggered = self.add_stagger_to_sequences(restored_sequences[3:])
        
        # Equalize sequence lengths
        template_staggered = self.equalize_sequence_lengths(template_staggered)
        input_staggered = self.equalize_sequence_lengths(input_staggered)
        
        staggered_restored_sequences = template_staggered + input_staggered

        # Step 6: Write output files
        # 6.1 Write MSA
        msa_output = f"{self.output_prefix}_alignment.fasta"
        with open(msa_output, 'w') as f:
            for i, seq in enumerate(template_staggered + input_staggered):
                group = "Template" if i < 3 else "Input"
                chain = chr(65 + (i % 3))  # A, B, or C
                f.write(f">{group}_{chain}\n{str(seq.seq)}\n")

        # 6.2 Write MODELLER-formatted alignment
        end_atom = self.extract_last_atom_serial_number(self.template_pdb)
        modeller_output = f"{self.output_prefix}_modeller.ali"
        self.write_modeller_formatted_output(staggered_restored_sequences, modeller_output, end_atom)

        return staggered_restored_sequences, msa_output, modeller_output

def align_sequences(input_fasta_path: str, template_fasta_path: str, output_prefix: str, template_pdb: str) -> Tuple[List[SeqRecord], str, str]:
    alignment = Alignment(input_fasta_path, template_fasta_path, output_prefix, template_pdb)
    return alignment.align_sequences()
