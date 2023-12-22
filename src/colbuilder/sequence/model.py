
class Model:
    """
    
    class represent the pdb file for the sequence alignment
    
    """
    def __init__(self,pdb_filename=None):
        self.pdb_filename=pdb_filename
        self.keys=['atom_id','atom_type','resname','modified_resname','chain_id','resid','x','y','z','element']
        self.pdb={ k:[] for k in self.keys }
    
    def read_pdb(self,pdb_filename=None):
        """
        
        read pdb and define as the model sequence
        
        """
        if pdb_filename==None: pdb_filename=self.pdb_filename

        with open(pdb_filename+'.pdb','r') as f:
            for l in f: 
                if l[0:4]=='ATOM':
                    line_split=[i for i in l.split(' ') if i!='']
                    self.pdb['atom_id'].append(line_split[1])
                    self.pdb['atom_type'].append(line_split[2])
                    self.pdb['resname'].append(line_split[3])
                    self.pdb['modified_resname'].append(line_split[3])
                    self.pdb['chain_id'].append(line_split[4][0])
                    
                    if len(line_split)==9:
                        self.pdb['resid'].append(line_split[4][1:])
                        self.pdb['x'].append(line_split[5])
                        self.pdb['y'].append(line_split[6])
                        self.pdb['z'].append(line_split[7])
                        self.pdb['element'].append(line_split[8])
                    elif len(line_split)==10:
                        self.pdb['resid'].append(line_split[5])
                        self.pdb['x'].append(line_split[6])
                        self.pdb['y'].append(line_split[7])
                        self.pdb['z'].append(line_split[8])
                        self.pdb['element'].append(line_split[9])
    
    def prepare_pdb(self,pdb=None):
        """
        
        prepare pdb-file for input in modeller
        
        """
        self.pdb['modified_resname']=[for resname in self.pdb['resname'] if resname=='HYP'