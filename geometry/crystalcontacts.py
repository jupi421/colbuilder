class CrystalContacts:
    """
    
    Reads contact information from crystal_contact file
    Writes updated contact information for chimera
    
    --
    
    input:  -f     crystal contacts file from chimera
   
    output: -o     optimized crystal contacts file for chimera      
        
    """
    def __init__(self,crystalcontact_file=None):
        self.crystalcontact_file=crystalcontact_file
        self.t_matrix={ }

    def read_contacts(self,crystalcontact_file=None):
        """
        
        Read crystal contacts information from chimera contact output-file
        
        """
        if crystalcontact_file==None: crystalcontact_file=self.crystalcontact_file
        return open(crystalcontact_file+'.txt','r').readlines()
    
    def read_t_matrix(self,crystalcontact_file=None,crystalcontacts=None):
        """
        
        Read transformation matrix T from contact file 
        
        """
        if crystalcontacts==None: crystalcontacts=self.read_contacts(crystalcontact_file)
        for idx in range(0,len(crystalcontacts),4):
            self.t_matrix[float(crystalcontacts[idx].split(' ')[1])]=[   
                float(crystalcontacts[idx+1].split(' ')[-1]),
                float(crystalcontacts[idx+2].split(' ')[-1]),
                float(crystalcontacts[idx+3].split(' ')[-1])]
        return self.t_matrix
    
    def write_contacts(self,system=None,crystalcontact_file=None):
        """
        
        Writes crystal contacts to txt file for chimera
        
        """
        if crystalcontact_file==None: crystalcontact_file=self.crystalcontact_file
        if system==None: 
            system=self.read_contacts(crystalcontact_file)
            with open(crystalcontact_file+'.txt','w') as f:
                for key_cc,val_cc in system:
                    f.write('Model '+str(float(key_cc))+'\n')
                    f.write('         1 0 0 %s\n' % (round(val_cc[0],3)))
                    f.write('         0 1 0 %s\n' % (round(val_cc[1],3)))
                    f.write('         0 0 1 %s\n' % (round(val_cc[2],3)))
            f.close()
        else:
            with open(crystalcontact_file+'.txt','w') as f:
                for model in system.get_keys():
                    f.write('Model '+str(float(model))+'\n')
                    f.write('         1 0 0 %s\n' % (system.get_model(model_id=model).model_t[0]))
                    f.write('         0 1 0 %s\n' % (system.get_model(model_id=model).model_t[1]))
                    f.write('         0 0 1 %s\n' % (system.get_model(model_id=model).model_t[2]))
            f.close()
        return

    def find_contact(self,model_id=None):
        """
        
        Finds the translation vector for one specific model-id
        
        """
        self.t_matrix=self.read_t_matrix(self.crystalcontact_file)
        return self.t_matrix[model_id]