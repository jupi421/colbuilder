class Crosslink:
    """
    
    Crosslink-class is a container for crosslink information for each model
    
    """
    def __init__(self,model_id=None,crosslink_id=None,resid=None,resname=None,chain=None):
        self.model_id=model_id
        self.id=crosslink_id
        self.resname={ }
        self.resid={ }
        self.chain={ }
    

