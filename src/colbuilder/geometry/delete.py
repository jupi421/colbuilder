import numpy as np

class Delete:
    """

    delete connected models with restrictions to obtain a less crosslinked system

    """
    def __init__(self,delete_ratio=None,system=None,fibril_length=None):
        self.delete=delete_ratio
        self.system=system
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
        self.z_min=self.system.get_model(model_id=0.0).cog - fibril_length/2
        self.z_max=self.system.get_model(model_id=0.0).cog + fibril_length/2

    def run_delete(self,delete_ratio=None,system=None):
        """
        
        set mutation of crystalcontacts according to deletion ratio
        
        """
        while self.system.count_states(state='delete') / ( self.system.count_states(state='none') + 
                                                    self.system.count_states(state='protect') + 
                                                    self.system.count_states(state='delete')  ) <  int(self.mutate_ratio) / 100:
            model=self.system.get_model(model_id=self.draw_model(system=self.system))
            self.delete_model(model=model)
        return self.system

    def delete_model(self,model=None):
        """

        change status of delete for crosslink and protect nearest neighbors

        """
        delete_id=self.draw_crosslink(model)
        if model.crosslink[delete_id].state!='delete' or model.crosslink[delete_id].state!='protect' or model.crosslink[delete_id].position[2]<self.z_min or model.crosslink[delete_id].position[2]>self.z_max: 
            model.crosslink[delete_id].state='delete'
            for cross in model.crosslink:
                if cross!=model.crosslink[delete_id] and np.linalg.norm(cross.position-model.crosslink[delete_id].position)<=10: cross.state='delete'
                elif cross!=model.crosslink[delete_id] and 11<=np.linalg.norm(cross.position-model.crosslink[delete_id].position)<=1000: cross.state='protect'
    
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
    
    def write_delete(self,system=None,mutation_file='delete'):
        """
        
        write deletions of system to file used as input for chimera
        
        """
        with open(mutation_file+'.txt','w') as f:
            for key in system.get_models():
                for cross in system.get_model(model_id=key).crosslink:
                    if cross.state=='delete': 
                        f.write(str(int(key))+'.caps.pdb '+str(cross.resname)+' '+str(cross.resid)+' '+str(cross.chain)+'\n' )
        f.close()