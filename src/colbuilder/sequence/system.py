import re

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
            'MSE':'M', 'HYP':'O', 'L4Y':'l', 'L5Y':'j'
            }
        self.pattern=re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
        self.crystal=''
    
    def read_pdb(self,pdb_filename=None):
        """
        
        read pdb and define as the model sequence
        
        """
        if pdb_filename==None: pdb_filename=self.pdb_filename
        cnt=0
        with open(pdb_filename+'.pdb','r') as f:
            for l in f:
                if l[0:4]=='ATOM':
                    line_list=[l[0:7],l[7:13],l[13:17],l[17:20],l[21:22],l[22:26],l[30:38],l[38:46],l[46:54],l[56:]]
                    self.atoms['atom_id'].append(line_list[1])
                    self.atoms['atom_type'].append(line_list[2])
                    self.atoms['atom_cnt'].append(str(cnt))
                    self.atoms['resname'].append(line_list[3])
                    self.atoms['chain_id'].append(line_list[4])
                    self.atoms['resid'].append(line_list[5])
                    self.atoms['x'].append(line_list[6])
                    self.atoms['y'].append(line_list[7])
                    self.atoms['z'].append(line_list[8])
                    self.atoms['element'].append(line_list[9].replace('\n',''))
                    cnt+=1
                elif l[0:5]=='CRYST':
                    self.crystal=l
    
    def pdb_to_fasta(self,atoms=None):
        """
        
        transform atoms of pdb file to fasta format
        
        """
        if atoms==None: atoms=self.atoms
        self.atoms['fasta']=[self.aminoacids[resname] for resname in atoms['resname']]