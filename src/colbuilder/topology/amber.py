import numpy as np

class Amber:
    """

    class to generate topology for the amber99 force field

    """
    def __init__(self,system=None,ff=None):
        self.system=system
        self.ff=ff
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
        self.is_write=('Protein_chain_A', '[ system ]\n')
    
    def merge_pdbs(self,connect_id=None):
        """
        
        merge pdb's according to connect_id in system
        
        """
        if self.system.get_model(model_id=connect_id).connect!=None:
            with open(str(int(connect_id))+'.merge.pdb','w') as f:
                for model in self.system.get_model(model_id=connect_id).connect:
                    pdb_model=open(str(self.system.get_model(model_id=connect_id).type)+'/'+str(int(model))+'.caps.pdb','r').readlines()
                    f.write("".join(i for i in pdb_model if i[0:6] in self.is_line))
                f.write("END")
            f.close()

    def write_itp(self,itp_file=None):
        """
        
        reads itp-file and cleans it
        
        """
        itp_model=open(str(itp_file),'r').readlines()
        write=0
        with open(str(itp_file),'w') as f:
            for i in itp_model:
                if i in self.is_write: 
                    f.write('[ moleculetype ]\n')
                    f.write(str(itp_file)+'  3')
                    write+=1
                if write==1:  f.write("".join(i for i in itp_model))
        f.close()
    
    def write_topology(self,system=None,topology_file=None):
        """
        
        writes a topology-file for amber99-ildnp-star
        
        """
        with open(topology_file,'w') as f:
            f.write('; Topology for Collagen Microfibril from Colbuilder 2.0')
            f.write('#include "./'+self.ff+'/forcefield.itp"\n')

            for model in system.get_keys():
                f.write('#include col_'+str(int(model))+'.itp\n')
            
            f.write('#include "./'+self.ff+'/ions.itp"\n')
            f.write('#include "./'+self.ff+'/tip3p.itp"\n')

            f.write('\n\n [ system ] \n ; name Collagen Microfibril from Colbuilder 2.0 \n\n [ molecules ] \n; name  number')

            for model in system.get_keys():
                f.write('col_'+str(int(model))+'   1\n')
    
    def write_gro(self,system=None,gro_file=None):
        """
        
        write a gro-file
        
        """
        gro=[]
        for model in system.keys():
            gro=open('col_'+str(int(model))+'.gro','r').readlines()
            for line in gro[2:-1]:
                with open('col_'+str(int(model))+'.gro','w') as f:
                    f.write(line)
                f.close()

