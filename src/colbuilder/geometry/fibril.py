import subprocess
import os
from colbuilder.geometry import model

class Fibril:
    """
    
    class to setup system from fibril
    
    """
    def __init__(self,pdb_file=None,system=None):
        self.pdb_file=pdb_file
        self.system=system
        self.fibril=open(pdb_file,'r').readlines()
        self.count=self.count_models(pdb_file=self.pdb_file)
        self.models={ k : {} for k in range(self.count) }

    def seperate_system(self,pdb_file=None):
        """
        
        separate system in models and write pdbs
        
        """
        if pdb_file==None: self.pdb_file=pdb_file ; self.fibril=open(pdb_file,'r').readlines()

        model=[] ; model_cnt=0
        for line in self.fibril:
            if line[0:3]=='TER' and last_line[21:22]=='C':
                self.write_model(model_cnt=model_cnt,lines=model) ; model=[] ; model_cnt+=1; continue
            model.append(line) ; last_line=line

    def write_model(self,model_cnt=None,lines=None):
        """
        
        writes model to pdb-file
        
        """
        with open(str(model_cnt)+'.caps.pdb','w') as f:
            for line in lines: f.write(line)
            f.write('END')
        f.close()

    def count_models(self,pdb_file=None):
        """
        
        counts number of models in fibril
        
        """
        cnt=0
        for line in self.fibril: 
            if 'TER' in line: cnt+=1
        return int(cnt/3)

    def build_system(self,system=None):
        """
        
        get a models from a pdb file and checks connection
        
        """
        if system==None: system=self.system
        for model_id in range(self.count):
            model_=model.Model(id=model_id,pdb_file=str(int(model_id))+'.caps')
            system.add_model(model=model_)
        return system

    def write_connect(self,system=None,connect_file=None):
        """
        
        writes system of model connections to file 
        
        """
        subprocess.run("rm -r N/ T/ D/ TD/ DT/ ",shell=True,
                       stdout=subprocess.DEVNULL,stderr=subprocess. DEVNULL)
        
        with open(connect_file+'.txt','w') as f:
            for model in system.get_models():
                if system.get_model(model_id=model).connect==None: continue
                elif len(system.get_model(model_id=model).connect)==1:
                    if not os.path.exists(os.getcwd()+'/N'): subprocess.run("mkdir N",shell=True)
                    f.write(str(int(model))+'.caps.pdb')
                    f.write(' ; N \n')
                    subprocess.run("mv "+str(int(model))+'.caps.pdb N/',shell=True,
                                   stdout=subprocess.DEVNULL,stderr=subprocess. DEVNULL)
                else:
                    if not os.path.exists(os.getcwd()+'/'+str( system.get_model(model_id=model).type)): 
                        subprocess.run("mkdir "+str( system.get_model(model_id=model).type),shell=True)

                    for connect in system.get_model(model_id=model).connect:
                        f.write(str(int(connect))+'.caps.pdb ')
                        subprocess.run("mv "+str(int(connect))+".caps.pdb "+str(system.get_model(model_id=model).type)+"/",shell=True)
                    f.write(' ; '+str(system.get_model(model_id=model).type)+'\n')

        f.close()