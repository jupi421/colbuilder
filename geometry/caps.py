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
    def __init__(self,system_size=None):
        self.system_size=system_size
        self.chains=['A','B','C']
        self.caps=['N','C']
        self.model={ k:[] for k in self.chains}

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
        return 

    def write_pymol(self,cap=None,chain_id=None):
        """
        
        Writes command to be used in pymol to add a cap

        """
        if cap=='N': index=0
        elif cap=='C': index=-1
        return 'resi '+str(self.model[chain_id][index])+' and chain '+str(chain_id)+' and name '+str(cap)
    
    def add_caps(self,pdb_id=None):
        """
        
        adds caps to both ends of each model
        
        """
        cmd.load(str(pdb_id)+'.pdb')
        for cap in self.caps:
            for chain in self.chains:
                line_cap=self.write_pymol(cap=cap,chain_id=chain)
                cmd.edit(line_cap)
                if cap=='N': editor.attach_amino_acid("pk1",'ace',ss=0)
                elif cap=='C': editor.attach_amino_acid("pk1",'nme',ss=0)
                break

