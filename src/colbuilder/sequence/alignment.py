# sequence_alignment.py

import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from io import StringIO

def align_sequences(input_fasta_path, template_fasta_path, output_prefix):
    alignment = Alignment(input_fasta_path, template_fasta_path, output_prefix)
    return alignment.align_sequences()

class Alignment:
    def __init__(self, input_fasta_path, template_fasta_path, output_prefix):
        self.input_fasta_path = input_fasta_path
        self.template_fasta_path = template_fasta_path
        self.output_prefix = output_prefix
        self.input_sequences = list(SeqIO.parse(input_fasta_path, "fasta"))
        self.template_sequences = list(SeqIO.parse(template_fasta_path, "fasta"))

    def add_stagger_to_sequences(self, sequences):
        for i, seq in enumerate(sequences):
            seq.seq = Seq('-' * i + str(seq.seq))
        return sequences

    def record_hydroxyproline_positions(self, sequence):
        return [pos for pos, char in enumerate(sequence) if char.upper() == "O"]

    def substitute_hydroxyproline(self, sequence, positions):
        sequence = list(sequence)
        for pos in positions:
            sequence[pos] = "P"
        return "".join(sequence)

    def create_position_mapping(self, original_seq, aligned_seq):
        map_original_to_aligned = {}
        aligned_index = 0
        for original_index in range(len(original_seq)):
            while aligned_index < len(aligned_seq) and aligned_seq[aligned_index] == "-":
                aligned_index += 1
            if aligned_index < len(aligned_seq):
                map_original_to_aligned[original_index] = aligned_index
                aligned_index += 1
        return map_original_to_aligned

    def restore_hydroxyproline_with_mapping(self, aligned_seq, original_positions, position_mapping):
        aligned_seq_list = list(aligned_seq)
        for original_pos in original_positions:
            if original_pos in position_mapping:
                new_pos = position_mapping[original_pos]
                aligned_seq_list[new_pos] = 'O'
        return "".join(aligned_seq_list)

    def get_muscle_version(self):
        result = subprocess.run(["muscle", "-version"], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            raise RuntimeError("MUSCLE is not installed or not found in PATH.")

    def align_sequences_with_muscle(self, input_path, output_path):
        muscle_version = self.get_muscle_version()
        if "3.8" in muscle_version:
            muscle_command = ["muscle", "-in", input_path, "-out", output_path]
        else:
            muscle_command = ["muscle", "-align", input_path, "-output", output_path]
        subprocess.run(muscle_command, check=True)
        print("Alignment completed.")

    def align_sequences(self):
        # Step 1: Align input sequences
        input_modified_sequences = []
        hydroxyproline_positions = {}
        
        for seq in self.input_sequences:
            positions = self.record_hydroxyproline_positions(seq.seq)
            hydroxyproline_positions[seq.id] = positions
            substituted_seq = self.substitute_hydroxyproline(seq.seq, positions)
            modified_seq = SeqRecord(Seq(substituted_seq), id=seq.id, description=seq.description)
            input_modified_sequences.append(modified_seq)

        temp_input_path = f"{self.output_prefix}_temp_input.fasta"
        temp_input_aligned_path = f"{self.output_prefix}_temp_input_aligned.afa"

        SeqIO.write(input_modified_sequences, temp_input_path, "fasta")
        self.align_sequences_with_muscle(temp_input_path, temp_input_aligned_path)
        
        input_aligned_sequences = list(SeqIO.parse(temp_input_aligned_path, "fasta"))

        # Step 2: Align the MSA from step 1 with the template
        combined_sequences = self.template_sequences + input_aligned_sequences
        temp_combined_path = f"{self.output_prefix}_temp_combined.fasta"
        temp_final_aligned_path = f"{self.output_prefix}_temp_final_aligned.afa"

        SeqIO.write(combined_sequences, temp_combined_path, "fasta")
        self.align_sequences_with_muscle(temp_combined_path, temp_final_aligned_path)

        final_aligned_sequences = list(SeqIO.parse(temp_final_aligned_path, "fasta"))

        # Step 3: Restore hydroxyprolines and add stagger
        restored_sequences = []
        for seq in final_aligned_sequences:
            if seq.id in hydroxyproline_positions:
                original_seq = [s for s in self.input_sequences if s.id == seq.id][0]
                position_mapping = self.create_position_mapping(str(original_seq.seq), str(seq.seq))
                restored_seq = self.restore_hydroxyproline_with_mapping(str(seq.seq), hydroxyproline_positions[seq.id], position_mapping)
                restored_sequences.append(SeqRecord(Seq(restored_seq), id=seq.id, description=seq.description))
            else:
                restored_sequences.append(seq)

        staggered_restored_sequences = self.add_stagger_to_sequences(restored_sequences)

        # Clean up temporary files
        for temp_file in [temp_input_path, temp_input_aligned_path, temp_combined_path, temp_final_aligned_path]:
            os.remove(temp_file)

        return staggered_restored_sequences
