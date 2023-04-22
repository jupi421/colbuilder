class Contacts:
    """
    
    Reads contact information from crystal_contact file
    Writes updated contact information for chimera
    
    --
    
    input:  -f     crystal contacts file from chimera
   
    output: -o     optimized crystal contacts file for chimera      
        
    """
    def __init__(self,contact_file=None):
        self.contact_file=contact_file
        self.t_matrix={ }

    
    def read_contacts(self,contact_file=None):
        """
        
        Read crystal contacts information from chimera contact output-file
        
        """
        if contact_file==None: contact_file=self.contact_file
        return open(contact_file+'.txt','r').readlines()
    
    def read_t_matrix(self,contact_file=None,crystal_contacts=None):
        """
        
        Read transformation matrix T from contact file 
        
        """
        if crystal_contacts==None: crystal_contacts=self.read_contacts(contact_file)
        for idx in range(0,len(crystal_contacts),4):
            self.t_matrix[float(crystal_contacts[idx].split(' ')[1])]=[   
                float(crystal_contacts[idx+1].split(' ')[-1]),
                float(crystal_contacts[idx+2].split(' ')[-1]),
                float(crystal_contacts[idx+3].split(' ')[-1])]
        return self.t_matrix
    
    def write_contacts(self,contact_file=None,crystal_contacts=None):
        """
        
        Writes crystal contacts to txt file for chimera
        
        """
        if contact_file==None: contact_file=self.contact_file
        if crystal_contacts==None: crystal_contacts=self.read_contacts(contact_file)
        with open(contact_file+'.txt','w') as f:
            for key_cc,val_cc in crystal_contacts:
                f.write(str(key_cc)+'\n')
                f.write('         1 0 0 %s\n' % (round(val_cc[0],3)))
                f.write('         1 0 0 %s\n' % (round(val_cc[1],3)))
                f.write('         1 0 0 %s\n\n' % (round(val_cc[2],3)))
        f.close()