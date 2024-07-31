# Adapted from https://github.com/kad-ecoli/pdb2fasta
import re
import sys

def pdb_to_fasta(pdb_file):
    aa3to1 = {
        "ALA": "A", "VAL": "V", "PHE": "F", "PRO": "P", "MET": "M",
        "ILE": "I", "LEU": "L", "ASP": "D", "GLU": "E", "LYS": "K",
        "ARG": "R", "SER": "S", "THR": "T", "TYR": "Y", "HIS": "H",
        "CYS": "C", "ASN": "N", "GLN": "Q", "TRP": "W", "GLY": "G",
        "MSE": "M", "HYP": "O", "L4Y": "-", "L5Y": "-"
    }
    
    ca_pattern = re.compile(
        "^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|"
        "^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA](HYP|L4Y|L5Y)\s([\s\w])|"
        "^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])"
    )
    
    chain_dict = {}
    chain_list = []
    
    with open(pdb_file, 'r') as fp:
        for line in fp:
            if line.startswith("ENDMDL"):
                break
            match_list = ca_pattern.findall(line)
            if match_list:
                resn = match_list[0][0] + match_list[0][2] + match_list[0][4]
                chain = match_list[0][1] + match_list[0][3] + match_list[0][5]
                if chain in chain_dict:
                    chain_dict[chain] += aa3to1[resn]
                else:
                    chain_dict[chain] = aa3to1[resn]
                    chain_list.append(chain)
    
    fasta_content = ""
    for chain in chain_list:
        fasta_content += f">{pdb_file}:{chain}\n{chain_dict[chain]}\n"
    
    return fasta_content

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python pdb2fasta.py <pdb_file>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    fasta_content = pdb_to_fasta(pdb_file)
    print(fasta_content)
