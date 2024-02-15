import re
import numpy as np

class System:
    """
    
    class represents methods to be performed on the pdb file to prepare for the sequence alignment step
    
    """
    def __init__(self,collagen_type=None,pdb_filename=None,register=None):
        self.pdb_filename=pdb_filename
        self.collagen_type=collagen_type
        self.register="".join([i for i in register])
        self.keys=['atom_id','atom_cnt','atom_type','resname','modified_resname','chain_id','resid','x','y','z','element','fasta']
        self.atoms={ k:[] for k in self.keys }
        self.aminoacids={ 
            'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
            'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
            'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
            'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
            'MSE':'M'
            }
        self.pattern=re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
    
    def read_pdb(self,pdb_filename=None):
        """
        
        read pdb and define as the model sequence
        
        """
        if pdb_filename==None: pdb_filename=self.pdb_filename
        cnt=0
        with open(pdb_filename+'.pdb','r') as f:
            for l in f: 
                if l[0:4]=='ATOM':
                    line_split=[i for i in l.split(' ') if i!='']
                    self.atoms['atom_id'].append(line_split[1])
                    self.atoms['atom_type'].append(line_split[2])
                    self.atoms['atom_cnt'].append(cnt)
                    self.atoms['resname'].append(line_split[3])
                    self.atoms['modified_resname'].append(line_split[3])
                    self.atoms['chain_id'].append(line_split[4][0])
                    cnt+=1
                    
                    if len(line_split)==9:
                        self.atoms['resid'].append(line_split[4][1:])
                        self.atoms['x'].append(line_split[5])
                        self.atoms['y'].append(line_split[6])
                        self.atoms['z'].append(line_split[7])
                        self.atoms['element'].append(line_split[8].replace('\n',''))
                    elif len(line_split)==10:
                        self.atoms['resid'].append(line_split[5])
                        self.atoms['x'].append(line_split[6])
                        self.atoms['y'].append(line_split[7])
                        self.atoms['z'].append(line_split[8])
                        self.atoms['element'].append(line_split[9].replace('\n',''))
    
    def prepare_pdb(self,atoms=None):
        """
        
        prepare pdb-file for input in modeller
        
        """
        self.atoms['modified_resname']=[ 'PRO' if resname=='HYP' else resname for resname in self.atoms['resname'] ]

    def pdb_to_fasta(self,atoms=None):
        """
        
        transform atoms of pdb file to fasta format
        
        """
        if atoms==None: atoms=self.atoms
        self.atoms['fasta']=[self.aminoacids[resname] for resname in atoms['modified_resname']]

    def write_pdb(self,pdb_filename=None,atoms=None):
        """
        
        write modified atoms of pdb file 
        
        """
        atom_id=self.atoms['atom_id']
        atom_type=self.atoms['atom_type']
        chain_id=self.atoms['chain_id']
        resname=self.atoms['resname']
        resid=self.atoms['resid']
        x=self.atoms['x']
        y=self.atoms['y']
        z=self.atoms['z']
        element=self.atoms['element']
        with open(pdb_filename+'_mod.pdb','w') as f:
            for i in range(len(atom_id)):
                f.write("%6s%5d%4s%3s%4d%8.3f%8.3f%8.3f%2s\n" % ('ATOM',int(atom_id[i]),atom_type[i],resname[i],chain_id[i],int(resid[i]),float(x[i]),float(y[i]),float(z[i]),element[i]))
            f.close()