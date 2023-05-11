import numpy as np

class Martini:
    """

    class to generate topology for the martini 3 force field

    """
    def __init__(self,system=None,ff=None):
        self.system=system
        self.ff=ff
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU' )
        self.is_chain=('A','B','C')

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
    
    def set_pdb(self,pdb_id=None):
        """
        
        prepare pdbs for Martinize2
        
        """
        pdb=open(str(int(pdb_id))+'.caps.pdb','r').readlines()

        first_cnt=int(pdb[1][22:26])-1
        cnt=0
        out=[]
        cnt_map=0
        map=[]
        for line in pdb:
            if line[0:3]=='END': break
            
            if line in self.is_line:
                if first_cnt<int(line[22:26]): 
                    first_cnt=int(line[22:26]); cnt+=1; cnt_map+=1
                if first_cnt>int(line[22:26]) and line[21:22] in self.is_chain: 
                    cnt=1; cnt_map+=1

                if cnt<10: out.append(line[:22]+'   '+str(int(cnt))+line[26:])
                if cnt<100: out.append(line[:22]+'  '+str(int(cnt))+line[26:])
                if cnt<1000: out.append(line[:22]+' '+str(int(cnt))+line[26:])
                if cnt<10000: out.append(line[:22]+''+str(int(cnt))+line[26:])

                if cnt_map<10: map.append(line[:22]+'   '+str(int(cnt))+line[26:])
                if cnt_map<100: map.append(line[:22]+'  '+str(int(cnt))+line[26:])
                if cnt_map<1000: map.append(line[:22]+' '+str(int(cnt))+line[26:])
                if cnt_map<10000: map.append(line[:22]+''+str(int(cnt))+line[26:])
        out.append(pdb[-1])

    def write_pdb(self,pdb=None,file=None):
        """
        
        writes pdb to file
        
        """
        with open(file,'w') as f:
            for l in pdb: f.write(l)
        f.close()

