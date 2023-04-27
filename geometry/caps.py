import os
from pymol import cmd, editor

class Caps:
    """

    Adding Caps to a single triple helix

    Original version implemented by Dr. Agnieszka Obarska-Kosinska from  
    
    Obarska-Kosinska A, Rennekamp B, Ünal A, Gräter F. 
    ColBuilder: A server to build collagen fibril models.
    Viophys J. 2021 Sep 7;120(17):3544-3549. 
    doi: 10.1016/j.bpj.2021.07.009. 
    Epub 2021 Jul 13. PMID: 34265261; PMCID: PMC8456305.


    """
    # TODO: Do not just copy Agnieszka awesomve version, but think to make it better
    def __init__(self,system=None,pdb=None):
        self.system=system
        self.pdb_file=pdb
        self.system_size=system.size_models
        self.chains=['A','B','C']
        self.caps=['N','C']
        self.chain_length={ k:0 for k in self.chains }
        self.model={ k:[] for k in self.chains}
        self.get_chain_length(system=system)

    def read_residues(self,pdb_id=None):
        """
        
        Reads pdb-file for each chain in the triple helix to obtain
        exact location of where to place the cap.

        --

        input   :   pdb file for single triple helix
        
        """
        with open(str(pdb_id)+'.pdb','r') as f:
            for l in f:
                if l[0:6] in ('ATOM  ', 'HETATM', 'ANISOU', 'TER   '):
                    if l[13:15]=='CA' and l[21]=='A': self.model['A'].append(int(l[22:26]))
                    elif l[13:15]=='CA' and l[21]=='B': self.model['B'].append(int(l[22:26]))
                    elif l[13:15]=='CA' and l[21]=='C': self.model['C'].append(int(l[22:26]))
        f.close()

    def write_pymol(self,cap=None,chain_id=None):
        """
        
        Writes command to be used in pymol to add a cap

        """
        if cap=='N': index=0
        elif cap=='C': index=-1
        return 'resi '+str(self.model[chain_id][index])+' and chain '+str(chain_id)+' and name '+str(cap)
    
    def get_chain_length(self,system=None):
        """
        
        Gets the chain length for the initial model 0.0

        """
        self.read_residues(pdb_id=system.crystal.pdb_file)
        for chain in self.model:
            self.chain_length[chain]=len(self.model[chain])

    def add_caps(self,pdb_id=None):
        """
        
        adds caps to both ends of each model
        
        """
        cmd.load(str(pdb_id)+'.pdb')
        for cap in self.caps:
            for chain in self.chains:
                line_cap=self.write_pymol(cap=cap,chain_id=chain)
                if int(line_cap.split(' ')[1])!=1 and cap=='N':
                    cmd.edit(line_cap)
                    print(line_cap.split(' ')[1])
                    editor.attach_amino_acid("pk1",'ace',ss=0)
                elif int(line_cap.split(' ')[1])!=int(self.chain_length[cap]) and cap=='C': 
                    cmd.edit(line_cap)
                    editor.attach_amino_acid("pk1",'nme',ss=0)
                break

