import modeller
from modeller.automodel import *
from modeller.automodel import AutoModel

class Modeller:
    """
    
    class prepares the aligned sequence to generate the triple helical structure
    
    """
    def __init__(self,system=None,sequence=None,file=None,fasta=None,ensemble=None):
        self.sequence={ k:[] for k in sequence.keys()}
        self.system=system
        self.fasta=fasta
        self.file=file
        self.ensemble=int(ensemble)
        self.score={}
        self.modeller_pdb=''
        modeller.log.minimal()
    
    def read_muscle(self,muscle_file=None):
        """

        read multiple-sequence-alignment-file from Muscle
        
        """
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
        with open(self.file+'.msa','w') as f:
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
        self.reorder(register=register)

        for k,v in self.sequence.items(): self.sequence[k]=[i.replace('\n','') for j in v for i in j if i !='']

    def check_alignment(self,align=None,alignment_file=None):
        """
        
        check the alignment from with the modeller env
        
        """
        align.append(file=alignment_file+'_modeller.ali',align_codes='all')
        align.write(file=alignment_file+'.pap',alignment_format='PAP')
        align.write(file=alignment_file+'.fasta',alignment_format='FASTA')
        align.check_sequence_structure(gapdist=3.9)

    def perform_alignment(self,align=None,alignment_file=None):
        """
        
        perform the alignment step with the modeller env
        
        """
        align.append(file=alignment_file+'_modeller.ali',align_codes='all')
        align.salign() 

        align.write(file=alignment_file+'_length.ali',alignment_format='PIR')
        align.write(file=alignment_file+'_length.pap',alignment_format='PAP')

    def run_modeller(self,system=None,alignment_file=None):
        """
        
        run the modeller software to perform the alignment
        
        """
        env_=modeller.Environ()
        env_.io.hetatm=True
        align_=modeller.Alignment(env=env_)

        self.write_alignment(alignment_file=alignment_file+'_modeller',system=system)
        self.check_alignment(align=align_,alignment_file=alignment_file)
        self.perform_alignment(align=align_,alignment_file=alignment_file)

        env_.libs.topology.read('${LIB}/top_heav.lib')
        env_.libs.parameters.read('${LIB}/par.lib')

        auto_model=AutoModel(env_,alnfile=alignment_file+'_length.ali',
                        knowns='template',sequence='target',
                        assess_methods=(assess.DOPE,assess.GA341))

        auto_model.starting_model=1
        auto_model.ending_model=self.ensemble

        auto_model.make()
        self.score={i['name']:i['GA341 score'][0] for i in auto_model.outputs}

        self.modeller_pdb=max(self.score,key=self.score.get)


    def write_alignment(self,system=None,alignment_file=None):
        """
        
        write input file for modeller software
        
        """
        cnt=0
        with open(alignment_file+'.ali','w') as f:
            f.write('>P1;template\n')
            f.write('structure:'+str(system.pdb_filename)+'_mod: '+str(system.collagen_type)+
                    ' :A:+'+str(int(system.atoms['atom_cnt'][-1])+3)+':C: : : :\n')
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