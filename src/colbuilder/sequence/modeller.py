import modeller 

class Modeller:
    """
    
    class prepares the aligned sequence to generate the triple helical structure
    
    """
    def __init__(self,system=None,sequence=None):
        self.sequence={ k:[] for k in sequence.keys()}
        self.system=system
    
    def read_muscle(self,muscle_file=None):
        """

        read multiple-sequence-alignment-file from Muscle
        
        """
        tmp=''
        with open(muscle_file,'r') as f:
            for line in f:
                if line[0]=='>': 
                    chain_key=line.split(':')[1].replace('\n','')
                    continue
                self.sequence[chain_key].append(line)
        f.close()
        return self.sequence
    
    def reorder(self,register=None):
        """
    
        reorders the alignment
    
        """
        with open('tmp.msa','w') as f:
            for reg_id in register:
                f.write('>tmp:'+str(reg_id)+'\n')
                for fasta_id in self.sequence[reg_id]: f.write(fasta_id.replace('\n',''))
                f.write('\n')
        f.close()

    def prepare_alignment(self,muscle_file=None,register=None):
        """
        
        prepare the alignment from faste as input for modeller
        
        """
        self.sequence=self.read_muscle(muscle_file=muscle_file)
        self.fasta=self.sequence
        self.reorder(register=register)

        for k,v in self.sequence.items(): self.sequence[k]=[i.replace('\n','') for j in v for i in j if i !='']

    def check_alignment(self,alignment_file=None):
        """
        
        check the alignment from with the modeller env
        
        """
        env_=modeller.Environ()
        env_.io.hetatm=True
        align_=modeller.Alignment(env=env_)
        align_.append(file=alignment_file,align_codes='all')
        align_.write(file=alignment_file.replace('ali','pap'),alignment_format='PAP')
        align_.write(file=alignment_file.replace('ali','fasta'),alignment_format='FASTA')
        align_.check_sequence_structure(gapdist=3.9)        


    def write_alignment(self,system=None,alignment_file=None):
        """
        
        write input file for modeller software
        
        """
        cnt=0
        with open(alignment_file,'w') as f:
            f.write('>P1;template\n')
            f.write('structure:'+str(system.pdb_filename)+'_mod: '+str(system.collagen_type)+
                    ':A:+'+str(system.atoms['atom_cnt'][-1]+3)+':C: : : :\n')
            for k in self.fasta:
                f.write("".join([v for v in self.fasta[k]]))
                cnt+=1
                if cnt<len(self.fasta): f.write('/')
                else: f.write('*')

            cnt=0
            f.write('\n')
            f.write('\n>P1;target\n')
            f.write('sequence:tmp'+str(system.register)+': : : : : : : :\n')
            for k in self.sequence:
                f.write("".join([v for v in self.sequence[k]]))
                cnt+=1
                if cnt<len(self.sequence): f.write('/')
                else: f.write('*')
        f.close()
            




