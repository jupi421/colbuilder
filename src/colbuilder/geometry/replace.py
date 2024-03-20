import numpy as np

class Replace:
    """

    remove crosslinks from ceratin models to reduce the overall number of crosslinks

    """
    def __init__(self,ratio_replace=None,system=None,fibril_length=None):
        self.ratio=ratio_replace
        self.system=system
        self.is_line=('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
        self.external={}
        self.z_min=self.system.get_model(model_id=0.0).get_cog() - fibril_length/2
        self.z_max=self.system.get_model(model_id=0.0).get_cog() + fibril_length/2

    def run_replace(self,ratio_replace=None,system=None):
        """
        
        select models with crosslinks to be replaced by lysines
        
        """
        if system==None: system=self.system

        while ( 

            system.count_states(state='replace') / ( 
            system.count_states(state='none') + 
            system.count_states(state='protect') + 
            system.count_states(state='replace')  
            ) < int(ratio_replace) / 100 
            
            ):

            model=system.get_model(model_id=self.draw_model(system=system))
            self.replace_model(model=model)
        return system

    def replace_model(self,model=None):
        """

        change status of replace for crosslink and protect nearest neighbors

        """
        id_=self.draw_crosslink(model)
        if ( 
            model.crosslink[id_].state!='replace' or model.crosslink[id_].state!='protect' or 
            model.crosslink[id_].position[2]<self.z_min or model.crosslink[id_].position[2]>self.z_max
            ): 
            model.crosslink[id_].state='replace'
            for cross in model.crosslink:
                if (
                    cross!=model.crosslink[id_] and 
                    np.linalg.norm(cross.position-model.crosslink[id_].position)<=10
                    ): 
                    cross.state='replace'
                elif (
                    cross!=model.crosslink[id_] and 
                    11<=np.linalg.norm(cross.position-model.crosslink[id_].position)<=1000
                    ): 
                    cross.state='protect'
    
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
    
    def write_replace(self,system=None,file='replace'):
        """
        
        write models with replacements to file, that is used as input for chimera
        
        """
        with open(file+'.txt','w') as f:
            for key in system.get_models():
                for cross in system.get_model(model_id=key).crosslink:
                    if cross.state=='replace': 
                        f.write(str(int(key))+'.caps.pdb '+str(cross.resname)+' '+str(cross.resid)+' '+str(cross.chain)+'\n' )
        f.close()
