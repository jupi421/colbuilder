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
    def __init__(self,id_model=None,t_matrix=None,s_matrix=None,contacts=None,fibril_id=None):
        self.id_model=id_model
        self.t_matrix=t_matrix
        self.s_matrix=s_matrix
        self.contacts=contacts
        self.fibril_id=fibril_id

    def add_contact_connect(self,contact_key=None,contact_connect=None):
        """
        
        Adds contacts information to each model
        
        """
        self.fibril_id=contact_key
        self.contacts=contact_connect[contact_key]
        
    
