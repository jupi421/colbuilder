
def build_model(model_id=None,connect_contacts=None,crystal=None,contacts=None):
    model={ k:None for k in ['id','shift_model','translate_model','contacts']}
    model['id']=model_id
    model['translate_model']=contacts.find_contact(model_id)
    model['shift_model']=crystal.get_s_matrix(t_matrix=model['translate_model'])
    model['contacts']=connect_contacts[model_id]
    return model

class Model:
    """
    
    Model is a container for each chain that stores all information:
    contact, crosslinks, ids, connections, transform_matrix, shift_matrix
    
    --
    id      - id of model
    s_node  - coordinates in unit-cell shift matrix space
    t_node  - cartesian coordinates in transformation matrix space
    edge    - ids of connected models

    """
    def __init__(self,id_model=None,crystal=None,contacts=None):
        self.id_model=id_model
        self.crystal=crystal
        self.contacts=contacts

    def run_build_model(self,model_id=None,connect_contacts=None,crystal=None,contacts=None):
        """
        
        Calls function to build model
        
        """
        if model_id==None: model_id=self.id_model
        if crystal==None: crystal=self.crystal
        if contacts==None: crystal=self.contacts
        return build_model(model_id=model_id,connect_contacts=connect_contacts,crystal=crystal,contacts=contacts)
    
    def run_build_contacts(self,connect_contacts=None,crystal=None,contacts=None):
        """
        
        Builds initial crystal contacts system by calling build model
        
        """
        if crystal==None: crystal=self.crystal
        if contacts==None: contacts=self.contacts
        models={ k:{ }  for k in connect_contacts}
        for model in connect_contacts:
            models[model]=self.run_build_model(model,connect_contacts,crystal,contacts)
        return models