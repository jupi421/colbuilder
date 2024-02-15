import subprocess

class Alignment:
    """
    
    class holds functions to perform the sequence alignment with fasta
    
    """
    def __init__(self,atoms=None):
        self.atoms=atoms 
        self.sequence={ k:[] for k in ['A','B','C']}
    
    def get_sequence(self,atoms=None):
        """
        
        get amino acid sequence from system
        
        """
        if atoms==None: atoms=self.atoms
        for atom_cnt in atoms['atom_cnt']:
            if atoms['atom_type'][atom_cnt]=='CA' and atoms['chain_id'][atom_cnt]=='A':
                self.sequence['A'].append(atoms['fasta'][atom_cnt])
            elif atoms['atom_type'][atom_cnt]=='CA' and atoms['chain_id'][atom_cnt]=='B':
                self.sequence['B'].append(atoms['fasta'][atom_cnt])
            elif atoms['atom_type'][atom_cnt]=='CA' and atoms['chain_id'][atom_cnt]=='C':
                self.sequence['C'].append(atoms['fasta'][atom_cnt])
    
    def write_sequence(self,sequence=None):
        """
        
        write sequence to file for MUSCLE
        
        """
        if sequence==None: sequence=self.sequence
        with open('tmp.fasta','w') as f:
            f.write('>tmp:A\n')
            for s in sequence['A']:
                f.write(str(s))
            f.write('\n>tmp:B\n')
            for s in sequence['B']:
                f.write(str(s))
            f.write('\n>tmp:C\n')
            for s in sequence['C']:
                f.write(str(s))
            f.write('\n')
        f.close()

    def align_sequence(self,atoms=None):
        """
        
        align sequence from pdb
        
        """
        self.get_sequence(atoms=atoms)
        self.write_sequence(sequence=self.sequence)

        subprocess.run('muscle -in tmp.fasta -out tmp.afa' ,shell=True,
                       stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)