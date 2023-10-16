from colbuilder.topology import crosslink
class Itp:
    """

    class to merge itps from the martini 3 force field and go-like potentials

    """
    def __init__(self,system=None,model_id=None):
        self.system=system
        self.molecule=self.allocate(model_id=model_id)
        self.mol_ends=self.allocate(model_id=model_id)
        self.atoms=self.allocate(model_id=model_id)
        self.posres=self.allocate(model_id=model_id)
        self.bonds=self.allocate(model_id=model_id)
        self.constraints=self.allocate(model_id=model_id)
        self.angles=self.allocate(model_id=model_id)
        self.exclusions=self.allocate(model_id=model_id)
        self.go_exclusions=self.allocate(model_id=model_id)
        self.dihedrals=self.allocate(model_id=model_id)
        self.virtual_sites=self.allocate(model_id=model_id)
        self.go_table=self.allocate(model_id=model_id)
        self.pairs=self.allocate(model_id=model_id)
        self.final_atoms=[]
        self.final_posres=[]
        self.final_bonds=[]
        self.final_flex_bonds=[]
        self.final_constraints=[]
        self.final_angles=[]
        self.final_dihedrals=[]
        self.final_exclusions=[]
        self.final_go_exclusions=[]
        self.final_virtual_sites=[]
        self.final_pairs=[]
        self.crosslink_bonds=[]
        self.vs_to_col={ }
        self.delta_merge=0
        self.no_line=('[','\n','#endif\n','#ifdef','#ifndef','#include',';[',';')

    def allocate(self,model_id=None):
        """
        
        allocate storage to be filled
        
        """
        size=len(self.system.get_model(model_id=model_id).connect)
        return [[] for s in range(size)]

    def read_model(self,model_id=None,system=None):
        """
        
        merge connected itps to single itp-file
        
        """
        cnt_con=0
        for connect_id in self.system.get_model(model_id=model_id).connect:
            self.read_itp(model_id=model_id,connect_id=connect_id,cnt_con=cnt_con)
            self.read_excl(model_id=model_id,connect_id=connect_id,cnt_con=cnt_con)
            self.read_table(model_id=model_id,connect_id=connect_id,cnt_con=cnt_con)
            cnt_con+=1

    def read_itp(self,model_id=None,connect_id=None,cnt_con=None):
        """
        
        read single itp file 
        
        """
        bonded_type=''
        with open('col_'+str(int(model_id))+'.'+str(int(connect_id))+'.itp','r') as f:
            for l in f:
                if l[0]==';': continue
                self.molecule[cnt_con].append(l)
                
                if l=='[ atoms ]\n': bonded_type='atoms'
                elif l=='[ position_restraints ]\n':
                    self.mol_ends[cnt_con]=int(self.molecule[cnt_con][len(self.molecule[cnt_con])-3].split(' ')[0])
                    bonded_type='posres'
                elif l=='[ bonds ]\n': bonded_type='bonds'
                elif l=='[ constraints ]\n': bonded_type='constraints'
                elif l=='[ virtual_sitesn ]\n': bonded_type='virtualsites'
                elif l=='[ angles ]\n': bonded_type='angles'
                elif l=='[ dihedrals ]\n': bonded_type='dihedrals'
                elif l=='[ exclusions ]\n': bonded_type='exclusions'
            
                if bonded_type=='atoms' and l.split(' ')[0] not in self.no_line:
                    self.atoms[cnt_con].append([k for k in l.split(' ') if all([k!='',k!='\n'])])
                elif bonded_type=='posres' and l.split(' ')[0] not in self.no_line:
                    self.posres[cnt_con].append([k for k in l.split(' ') if all([k!='',k!='\n'])])
                elif bonded_type=='bonds' and l.split(' ')[0] not in self.no_line:
                    self.bonds[cnt_con].append([k for k in l.split(' ') if all([k!='',k!='\n'])])
                elif bonded_type=='constraints' and l.split(' ')[0] not in self.no_line:
                    self.constraints[cnt_con].append([k for k in l.split(' ') if all([k!='',k!='\n'])])
                elif bonded_type=='virtualsites' and l.split(' ')[0] not in self.no_line:
                    self.virtual_sites[cnt_con].append([k for k in l.split(' ') if all([k!='',k!='\n'])])
                elif bonded_type=='angles' and l.split(' ')[0] not in self.no_line:
                    self.angles[cnt_con].append([k for k in l.split(' ') if all([k!='',k!='\n'])])
                elif bonded_type=='exclusions' and l.split(' ')[0] not in self.no_line:
                    self.exclusions[cnt_con].append([k for k in l.split(' ') if all([k!='',k!='\n'])])
                elif bonded_type=='dihedrals' and l.split(' ')[0] not in self.no_line:
                    self.dihedrals[cnt_con].append([k for k in l.split(' ') if all([k!='',k!='\n'])])
        f.close()
    
    def read_excl(self,model_id=None,connect_id=None,cnt_con=None):
        """
        
        read go-excl file 
        
        """
        with open('col_'+str(int(model_id))+'.'+str(int(connect_id))+'_go-excl.itp','r') as f:
            for l in f:
                if l.split(' ')[0] not in self.no_line:
                    self.go_exclusions[cnt_con].append([k for k in l.split(' ') if all([k!='',k!='\n',k.strip()!=';'])])
        f.close()

    def read_table(self,model_id=None,connect_id=None,cnt_con=None):
        """
        
        read go-table file 
        
        """
        with open('col_'+str(int(model_id))+'.'+str(int(connect_id))+'_go-table.itp','r') as f:
            for l in f:
                if l.split(' ')[0] not in self.no_line:
                    self.go_table[cnt_con].append([k for k in l.split(' ') if all([k!='',k!='\n',k.strip()!=';'])])
        f.close()

    def go_to_pairs(self,model_id=None):
        """
        
        prepare pairs from go-like potentials
        
        """
        for cnt_con in range(len(self.system.get_model(model_id=model_id).connect)):
            self.match_vs_to_pairs(cnt_con=cnt_con)
            self.get_pairs(cnt_con=cnt_con)
    
    def match_vs_to_pairs(self,cnt_con=None):
        """
        
        match virtual sites from Go-model to pairs description
        
        """
        for a in range(len(self.atoms[cnt_con])):
            if self.atoms[cnt_con][a][1][0:3]=='col':
                self.vs_to_col[self.atoms[cnt_con][a][1]]=self.atoms[cnt_con][a][0]
                self.atoms[cnt_con][a][1]='col'
    
    def get_pairs(self,cnt_con=None):
        """
        
        get pairs to define the topology for the crosslinked molecules 
        
        """
        for t in range(len(self.go_table[cnt_con])):
            self.pairs[cnt_con].append([self.vs_to_col[self.go_table[cnt_con][t][0]],
            self.vs_to_col[self.go_table[cnt_con][t][1]],self.go_table[cnt_con][t][2],
            self.go_table[cnt_con][t][3],self.go_table[cnt_con][t][4],self.go_table[cnt_con][t][5],
            self.go_table[cnt_con][t][6],self.go_table[cnt_con][t][7],self.go_table[cnt_con][t][0],
            self.go_table[cnt_con][t][1]])

    def merge_topology(self,cnt_con=None):
        """
        
        set the final merged topology file
        
        """
        if cnt_con!=0: self.delta_merge+=self.mol_ends[cnt_con-1]

        tmp=[[int(k[0])+self.delta_merge,k[1],k[2],k[3],k[4],int(k[5])+self.delta_merge,k[6]] for k in self.atoms[cnt_con]]
        for k in tmp: self.final_atoms.append(k)

        tmp=[[int(k[0])+self.delta_merge,k[1],k[2],k[3],k[4]] for k in self.posres[cnt_con]]
        for k in tmp: self.final_posres.append(k)

        tmp=[[int(k[0])+self.delta_merge,int(k[1])+self.delta_merge,k[2],k[3],k[4]] for k in self.bonds[cnt_con]]
        for k in tmp:
            if int(k[-1])==1000000: self.final_flex_bonds.append(k)
            else: self.final_bonds.append(k)
        
        tmp=[[int(k[0])+self.delta_merge,int(k[1])+self.delta_merge,int(k[2])+self.delta_merge,k[3],k[4],k[5]] for k in self.angles[cnt_con]]
        for k in tmp: self.final_angles.append(k)
    
        if len(self.dihedrals[cnt_con][1])==8:
            tmp=[[int(k[0])+self.delta_merge,int(k[1])+self.delta_merge,int(k[2])+self.delta_merge,int(k[3])+self.delta_merge,k[4],k[5],k[6],k[7]] for k in self.dihedrals[cnt_con] if len(k)==8]
        elif len(self.dihedrals[cnt_con][1])==7:
            tmp=[[int(k[0])+self.delta_merge,int(k[1])+self.delta_merge,int(k[2])+self.delta_merge,int(k[3])+self.delta_merge,k[4],k[5],k[6]] for k in self.dihedrals[cnt_con] if len(k)==7]
        elif len(self.dihedrals[cnt_con][1])==6:
            tmp=[[int(k[0])+self.delta_merge,int(k[1])+self.delta_merge,int(k[2])+self.delta_merge,int(k[3])+self.delta_merge,k[4],k[5]] for k in self.dihedrals[cnt_con] if len(k)==6]
        for k in tmp: self.final_dihedrals.append(k)

        tmp=[[int(k[0])+self.delta_merge,int(k[1])+self.delta_merge,k[2],k[3]] for k in self.constraints[cnt_con]]
        for k in tmp: self.final_constraints.append(k)
    
        tmp=[[str(int(k[0])+self.delta_merge),str(int(1)),str(int(k[2])+self.delta_merge)+'\n'] for k in self.virtual_sites[cnt_con]]
        for k in tmp: self.final_virtual_sites.append(k)

        tmp=[[int(k[0])+self.delta_merge,int(k[1])+self.delta_merge,k[2],int(k[3])+self.delta_merge,int(k[4])+self.delta_merge] for k in self.go_exclusions[cnt_con]]
        for k in tmp: self.final_go_exclusions.append(k)
    
        tmp=[]
        for k in self.exclusions[cnt_con]: tmp.append([str(int(j)+self.delta_merge) for j in k if j!='']) ; tmp[-1][-1]+='\n'
        for k in tmp: self.final_exclusions.append(k)
    
        tmp=[[int(k[0])+self.delta_merge,int(k[1])+self.delta_merge,k[2],format(float(k[3]),'.10f'),format(float(k[4]),'.10f'),';',k[-2],k[-1]+'\n'] for k in self.pairs[cnt_con]]
        for k in tmp: self.final_pairs.append(k)

    def make_topology(self,model_id=None,cnt_model=None):
        """
        
        make topology for merged itps and crosslinks
        
        """
        self.crosslink_bonds=crosslink.Crosslink(model_id=model_id).set_crosslink_bonds(model_id=model_id)
        for cnt_con in range(len(self.system.get_model(model_id=model_id).connect)):
            self.merge_topology(cnt_con=cnt_con)
        self.write_topology(cnt_model=cnt_model)
        self.write_excl(cnt_model=cnt_model)

    def write_topology(self,cnt_model=None):
        """
        
        write merged topology to itp file
        
        """
        with open('col_'+str(int(cnt_model))+'.itp','w') as f:
            f.write('; Merging of topologies for models due to system\n')
            f.write('[ moleculetype ]\n')
            f.write('col_'+str(cnt_model)+' 1\n')
            f.write('\n\n[ atoms ]\n')
            for a in self.final_atoms: f.write('{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}\n'.format(a[0],a[1],a[2],a[3],a[4],a[5],a[6]))

            f.write('\n[ position_restraints ]\n')
            f.write('#ifdef POSRES\n')
            f.write(" ".join(str(i) for a in self.final_posres for i in a))
    
            f.write('#endif\n')
            f.write('\n[ bonds ]\n')
            f.write(" ".join(str(i) for b in self.final_bonds for i in b))
    
            f.write('#ifdef FLEXIBLE\n; side chain flexible\n')
            f.write(" ".join(str(i) for sc in self.final_flex_bonds for i in sc))
            f.write('#endif\n')
    
            f.write('; crosslink bonds \n')
            f.write(" ".join(str(i) for cb in self.crosslink_bonds for i in cb))
    
            f.write('\n[ pairs ]\n')
            f.write(" ".join(str(i) for p in self.final_pairs for i in p))
    
            f.write('\n[ constraints ]\n')
            f.write('#ifndef FLEXIBLE\n')
            f.write(" ".join(str(i) for c in self.final_constraints for i in c))
            f.write('#endif\n')
    
            f.write('\n[ virtual_sitesn ]\n')
            f.write(" ".join(str(i) for v in self.final_virtual_sites for i in v))
    
            f.write('\n[ angles ]\n')
            f.write(" ".join(str(i) for a in self.final_angles for i in a))
    
            f.write('\n[ dihedrals ]\n')
            f.write(" ".join(str(i) for d in self.final_dihedrals for i in d))
    
            f.write('\n[ exclusions ]\n')
            f.write(" ".join(str(i) for ex in self.final_exclusions for i in ex))
        f.close()

    def write_excl(self,cnt_model=None):
        """
        
        write merged exclusion to itp file
        
        """
        with open('col_'+str(cnt_model)+'_go-excl.itp','w') as f:
            f.write(';[ exclusions ]\n')
            for ex in self.final_go_exclusions:
                for k in ex[:-1]: f.write(str(k)+' ')
                f.write('\n')
        f.close()