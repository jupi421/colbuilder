import numpy as np

class Mutate:
    """

    mutate connected models with restrictions to obtain mutated crosslinked system

    """
    def __init__(self,mutate_ratio=None,system=None,fibril_length=None):
        self.mutate_ratio=mutate_ratio
        self.system=system
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
        self.z_min=self.system.get_model(model_id=0.0).cog - fibril_length/2
        self.z_max=self.system.get_model(model_id=0.0).cog + fibril_length/2

    def run_mutate(self,mutate_ratio=None,system=None):
        """
        
        set mutation of crystalcontacts according to mutate ratio
        
        """
        while self.system.count_states(state='mut') / ( self.system.count_states(state='no') + 
                                                    self.system.count_states(state='prot') + 
                                                    self.system.count_states(state='mut')  ) <  int(self.mutate_ratio) / 100:
            model=self.system.get_model(model_id=self.draw_model(system=self.system))
            self.mutate_model(model=model)
        return self.system

    def mutate_model(self,model=None):
        """

        change status of mutate for crosslink and protects nearest neighbor

        """
        mut_id=self.draw_crosslink(model)
        if model.crosslink[mut_id].state!='mut' or model.crosslink[mut_id].state!='prot' or model.crosslink[mut_id].position[2]<self.z_min or model.crosslink[mut_id].position[2]>self.z_max: 
            model.crosslink[mut_id].state='mut'
            for cross in model.crosslink:
                if cross!=model.crosslink[mut_id] and np.linalg.norm(cross.position-model.crosslink[mut_id].position)<=10: cross.state='mut'
                elif cross!=model.crosslink[mut_id] and 11<=np.linalg.norm(cross.position-model.crosslink[mut_id].position)<=1000: cross.state='prot'
    
    def draw_model(self,system):
        """
        
        draw a random model from system
        
        """
        return np.random.choice(system.get_models())
    
    def draw_crosslink(self,model):
        """
        
        draw a random crosslink from model
        
        """
        return np.random.randint(0,len(model.crosslink))
    
    def write_mutate(self,system=None,mutation_file='mutation'):
        """
        
        write mutations of system to file used as input for chimera
        
        """
        with open(mutation_file+'.txt','w') as f:
            for key in system.get_models():
                for cross in system.get_model(model_id=key).crosslink:
                    if cross.state=='mut': f.write(str(int(key))+'.caps.pdb '+str(cross.resname)+' '+str(cross.resid)+' '+str(cross.chain)+'\n' )
        f.close()