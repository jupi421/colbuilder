class Martini:
    """

    class to generate topology for the martini 3 force field

    """
    def __init__(self,system=None,force_field=None):
        self.system=system
        self.ff=force_field
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU')
        self.is_chain=('A','B','C')

    def merge_pdbs(self,model_id=None,cnt_model=None):
        """
        
        merge pdb's according to connect_id in system
        
        """
        if self.system.get_model(model_id=model_id).connect!=None:
            with open(str(int(cnt_model))+'.merge.pdb','w') as f:
                for connect_id in self.system.get_model(model_id=model_id).connect:
                    pdb_model=open(str(int(model_id))+'.'+str(int(connect_id))+'.CG.pdb','r').readlines()
                    f.write("".join(i for i in pdb_model if i[0:6] in self.is_line))
                f.write("END")
            f.close()
    
    def read_pdb(self,pdb_id=None):
        """
        
        read pdb for Martinize2
        
        """
        pdb=[]
        with open(self.system.get_model(model_id=pdb_id).type+'/'+str(int(pdb_id))+'.caps.pdb','r') as file:
            for line in file:
                if line[0:6] in self.is_line:
                    pdb.append(line)
        return pdb
    
    def set_pdb(self,pdb=None):
        """
        
        prepare pdbs for Martinize2
        
        """
        first_cnt=int(pdb[1][22:26])-1
        cnt,cnt_map=0,0
        order,map=[],[]
        chain_store='A'
        for line in pdb:

            if line[21:22]!=chain_store:
                order.append('TER\n')
                if chain_store=='A': chain_store='B'
                elif chain_store=='B': chain_store='C'

            if first_cnt<int(line[22:26]): 
                first_cnt=int(line[22:26]); cnt+=1; cnt_map+=1
            if first_cnt>int(line[22:26]) and str(line[21:22]) in self.is_chain: 
                first_cnt=int(line[22:26]); cnt=1; cnt_map+=1

            if cnt<10: order.append(line[:22]+'   '+str(int(cnt))+line[26:])
            elif 10<=cnt<100: order.append(line[:22]+'  '+str(int(cnt))+line[26:])
            elif 100<=cnt<1000: order.append(line[:22]+' '+str(int(cnt))+line[26:])
            elif 1000<=cnt<10000: order.append(line[:22]+str(int(cnt))+line[26:])

            if cnt_map<10: map.append(line[:22]+'   '+str(int(cnt_map))+line[26:])
            elif 10<=cnt_map<100: map.append(line[:22]+'  '+str(int(cnt_map))+line[26:])
            elif 100<=cnt_map<1000: map.append(line[:22]+' '+str(int(cnt_map))+line[26:])
            elif 1000<=cnt_map<10000: map.append(line[:22]+str(int(cnt_map))+line[26:])

        return order,map

    def translate_pdb(self,pdb=None):
        """
        
        translate pdbs for Martinize2
        
        """
        return [line[:46]+'{:.3f}'.format(round(float(line[46:54]),3))+line[54:] for line in pdb if line[0:6] in self.is_line]

    def cap_pdb(self,pdb=None):
        """
        
        cap pdbs according to connect_id in system
        
        """
        cter,nter='NME','ACE'
        chain_length=self.get_chain_length(pdb)
        for line_it in range(len(pdb)):
            if pdb[line_it][17:20]=='ALA':
                if pdb[line_it][21:26]=='A'+chain_length['A']:
                    pdb[line_it]=pdb[line_it][0:17]+'CLA '+pdb[line_it][21:]
                elif pdb[line_it][21:26]=='B'+chain_length['B']:
                    pdb[line_it]=pdb[line_it][0:17]+'CLA '+pdb[line_it][21:]
                elif pdb[line_it][21:26]=='C'+chain_length['C']:
                    pdb[line_it]=pdb[line_it][0:17]+'CLA '+pdb[line_it][21:]

        if pdb[2][17:20]=='GLN': nter='N-ter'
        if pdb[-2][17:20]=='CLA': cter='CLA'
        return pdb,cter,nter

    def get_chain_length(self,pdb=None):
        """
        
        Get length for each chain of the triple helix
        
        """
        chain_length={ key:'' for key in ['A','B','C']}
        for line_it in range(len(pdb)):
            if pdb[line_it][21:22]=='A' and pdb[line_it+1][21:22]=='B':
                chain_length['A']=pdb[line_it][22:26]
            if pdb[line_it][21:22]=='B' and pdb[line_it+1][21:22]=='C':
                chain_length['B']=pdb[line_it][22:26]
            if pdb[line_it][21:22]=='C' and line_it==int(len(pdb)-1):
                chain_length['C']=pdb[line_it][22:26]
        return chain_length

    def write_pdb(self,pdb=None,file=None):
        """
        
        writes pdb to file
        
        """
        with open(file,'w') as f:
            for l in pdb: 
                if l[0:3]!='END': f.write(l)
            f.write('END')
        f.close()
    
    def get_system_pdb(self,size=None):
        """
        
        write system pdb 
        
        """
        pdb=[]
        for cnt_model in range(size):
            with open(str(int(cnt_model))+'.merge.pdb','r') as f:
                for l in f: pdb.append(l)
            f.close()
        return pdb  
      
    def write_system_topology(self,size=None):
        """
        
        write final topology file 
        
        """
        with open('system.top','w') as f:
            f.write('; This is the topology for the collagen microibril\n')
            f.write('#define GO_VIRT\n')
            f.write('#include "martini_v3.0.0.itp"\n')
            f.write('#include "go-sites.itp"\n\n')

            for m in range(size):
                f.write('#include "col_'+str(m)+'.itp"\n')
                f.write('#include "col_'+str(m)+'_go-excl.itp"\n')
            f.write('\n#include "martini_v3.0.0_solvents_v1.itp"\n')
            f.write('#include "martini_v3.0.0_ions_v1.itp"\n')

            f.write('\n[ system ]\n')
            f.write('Collagen, Martini 3 and Go-Potentials \n')
            f.write('\n[ molecules ]\n')
            for t in range(size):
                f.write('col_'+str(t)+'     1\n')
        f.close()

        self.write_go_topology(name_type='sites.itp')

    def write_go_topology(self,name_type=None):
        """
        
        write topology for go-like potentials
        
        """
        open('go-'+name_type,'w').close()
        with open('go-'+name_type,'w') as f:
            f.write('#include "col_go-'+name_type+'"\n')
        f.close()

        with open('col_go-'+name_type,'w') as f:
            f.write('[ atomtypes ]\n')
            f.write('; protein BB virtual particle\n')
            f.write('col 0.0 0.000 A 0.0 0.0 \n')
        f.close()
