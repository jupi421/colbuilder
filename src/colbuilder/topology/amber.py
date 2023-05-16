import os
import subprocess

class Amber:
    """

    class to generate topology for the amber99 force field

    """
    def __init__(self,system=None,ff=None):
        self.system=system
        self.ff=ff+'.ff'
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
    
    def merge_pdbs(self,connect_id=None):
        """
        
        merge pdb's according to connect_id in system
        
        """
        if self.system.get_model(model_id=connect_id).connect==None: return None
        if len(self.system.get_model(model_id=connect_id).connect)>1:
            with open(str(self.system.get_model(model_id=connect_id).type)+'/'+str(int(connect_id))+'.merge.pdb','w') as f:
                for model in self.system.get_model(model_id=connect_id).connect:
                    pdb_model=open(str(self.system.get_model(model_id=connect_id).type)+'/'+str(int(model))+'.caps.pdb','r').readlines()
                    f.write("".join(i for i in pdb_model if i[0:6] in self.is_line))
                f.write("END")
            f.close()
            return str(self.system.get_model(model_id=connect_id).type)

        elif len(self.system.get_model(model_id=connect_id).connect)==1:
            if os.path.exists("N/"+str(int(connect_id))+".caps.pdb"):
                subprocess.run("mv N/"+str(int(connect_id))+".caps.pdb N/"+str(int(connect_id))+".merge.pdb",
                           shell=True)
            return "N"


    def write_itp(self,itp_file=None):
        """
        
        reads itp-file and cleans it
        
        """
        itp_model=open(str(itp_file),'r').readlines()
        subprocess.run("rm "+str(itp_file),shell=True)

        write=False
        with open(str(itp_file).replace("top","itp"),'w') as f:
            for line in itp_model:
                if 'Include water topology' in line: break
                if write==True:  f.write(line)
                elif 'Protein_chain_A' in line: 
                    f.write('[ moleculetype ]\n')
                    f.write(str(itp_file).replace(".top","")+'  3\n')
                    write=True
        f.close()
    
    def write_topology(self,system=None,topology_file=None):
        """
        
        writes a topology-file for amber99-ildnp-star
        
        """
        with open(topology_file,'w') as f:
            f.write('; Topology for Collagen Microfibril from Colbuilder 2.0\n')
            f.write('#include "./'+self.ff+'/forcefield.itp"\n')

            for model in system.get_models():
                if os.path.exists("col_"+str(int(model))+".itp"):
                    f.write('#include "col_'+str(int(model))+'.itp"\n')
            
            f.write('#include "./'+self.ff+'/ions.itp"\n')
            f.write('#include "./'+self.ff+'/tip3p.itp"\n')

            f.write('\n\n[ system ]\n ;name\nCollagen Microfibril in Water\n\n[ molecules ]\n;name  number\n')

            for model in system.get_models():
                if os.path.exists("col_"+str(int(model))+".itp"):
                    f.write('col_'+str(int(model))+'   1\n')
    
    def write_gro(self,system=None,gro_file=None):
        """
        
        write a gro-file
        
        """
        gro=[]
        with open(gro_file,'w') as f:
            f.write("GROMACS GRO-FILE")
            for model in system.get_models():
                if os.path.exists("col_"+str(int(model))+".gro"):
                    gro=open('col_'+str(int(model))+'.gro','r').readlines()
                    for line in gro[2:-1]: f.write(line)
                    subprocess.run("rm col_"+str(int(model))+".gro",shell=True)
            f.write(gro[-1])
        f.close()

        len_gro=len(open(gro_file,'r').readlines())-1
        subprocess.run("sed -i '2i "+str(len_gro)+"' "+str(gro_file),shell=True)