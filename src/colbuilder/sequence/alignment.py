# sequence_alignment.py

import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO

class Alignment:
    def __init__(self, fasta_content, file_prefix):
        self.fasta_content = fasta_content
        self.file_prefix = file_prefix
        self.sequences = list(SeqIO.parse(StringIO(fasta_content), "fasta"))

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
            while aligned_seq[aligned_index] == "-":
                aligned_index += 1
            map_original_to_aligned[original_index] = aligned_index
            aligned_index += 1
        return map_original_to_aligned

    def restore_hydroxyproline_with_mapping(self, aligned_seq, original_positions, position_mapping):
        aligned_seq_list = list(aligned_seq)
        for original_pos in original_positions:
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

    def reorder_sequences_to_original(self, aligned_sequences, original_sequences):
        ordered_sequences = []
        for orig_seq in original_sequences:
            for aligned_seq in aligned_sequences:
                if orig_seq.id == aligned_seq.id:
                    ordered_sequences.append(aligned_seq)
                    break
        return ordered_sequences

    def align_sequences(self):
        modified_sequences = []
        for seq in self.sequences:
            positions = self.record_hydroxyproline_positions(seq.seq)
            substituted_seq = self.substitute_hydroxyproline(seq.seq, positions)
            modified_seq = SeqRecord(Seq(substituted_seq), id=seq.id, description=seq.description)
            modified_sequences.append(modified_seq)

        temp_input_path = f"{self.file_prefix}_temp_input.fasta"
        temp_output_path = f"{self.file_prefix}_temp_aligned.afa"

        SeqIO.write(modified_sequences, temp_input_path, "fasta")
        self.align_sequences_with_muscle(temp_input_path, temp_output_path)
        
        aligned_sequences = list(SeqIO.parse(temp_output_path, "fasta"))
        ordered_aligned_sequences = self.reorder_sequences_to_original(aligned_sequences, self.sequences)

        restored_sequences = []
        for seq in ordered_aligned_sequences:
            original_seq = [s for s in self.sequences if s.id == seq.id][0]
            positions = self.record_hydroxyproline_positions(original_seq.seq)
            position_mapping = self.create_position_mapping(str(original_seq.seq), str(seq.seq))
            restored_seq = self.restore_hydroxyproline_with_mapping(str(seq.seq), positions, position_mapping)
            restored_sequences.append(SeqRecord(Seq(restored_seq), id=seq.id, description=seq.description))

        staggered_restored_sequences = self.add_stagger_to_sequences(restored_sequences)

        # Clean up temporary files
        os.remove(temp_input_path)
        os.remove(temp_output_path)

        return staggered_restored_sequences

def align_sequences(fasta_content, file_prefix):
    alignment = Alignment(fasta_content, file_prefix)
    return alignment.align_sequences()
    
