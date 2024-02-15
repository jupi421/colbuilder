import subprocess

class Alignment:
    """
    
    class holds functions to perform the sequence alignment with fasta
    
    """
    def __init__(self,atoms=None,file=None):
        self.atoms=atoms 
        self.sequence={ k:[] for k in ['A','B','C']}
        self.file=file
    
    def get_sequence(self,atoms=None):
        """
        
        get amino acid sequence from system
        
        """
        if atoms==None: atoms=self.atoms
        for atom_cnt in atoms['atom_cnt']:
            if atoms['atom_type'][int(atom_cnt)].replace(' ','')=='CA' and atoms['chain_id'][int(atom_cnt)]=='A':
                self.sequence['A'].append(atoms['fasta'][int(atom_cnt)])
            elif atoms['atom_type'][int(atom_cnt)].replace(' ','')=='CA' and atoms['chain_id'][int(atom_cnt)]=='B':
                self.sequence['B'].append(atoms['fasta'][int(atom_cnt)])
            elif atoms['atom_type'][int(atom_cnt)].replace(' ','')=='CA' and atoms['chain_id'][int(atom_cnt)]=='C':
                self.sequence['C'].append(atoms['fasta'][int(atom_cnt)])
    
    def write_sequence(self,sequence=None):
        """
        
        write sequence to file for MUSCLE
        
        """
        if sequence==None: sequence=self.sequence
        with open(self.file+'.fasta','w') as f:
            f.write('>'+self.file+':A\n')
            for s in sequence['A']:
                f.write(str(s))
            f.write('\n>'+self.file+':B\n')
            for s in sequence['B']:
                f.write(str(s))
            f.write('\n>'+self.file+':C\n')
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

        subprocess.run('muscle -in '+self.file+'.fasta -out '+self.file+'.afa' ,
                       shell=True)#,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    
        self.fasta=self.sequence