class Model:
    """
    
    Class for each model in the system. Each model contains all information 
    associdated to id:

    --

    input   :   

                - model-id
                - translation vector
                - shift-matrix
                - neighbor contacts
                - fibril-id
                - ...

    output  :   class : model

    """
    def __init__(self,model_id=None,model_t=None,model_s=None,model_connect=None,model_connect_id=None):
        self.model_id=model_id
        self.model_t=model_t
        self.model_s=model_s
        self.model_connect=model_connect
        self.model_connect_id=model_connect_id

    def add_model_connect(self,model_connect_id=None,model_connect=None):
        """
        
        Adds information about connection of each model
        
        """
        self.model_connect_id=model_connect_id
        self.model_connect=model_connect
    
