import subprocess
import os
import numpy as np

class Chimera:
    """
    
    Generate collagen microfibril with the "crystal contacts" command from UCSF 
    Chimera. Chimera needs to be installed beforehand from:
    
    https://www.cgl.ucsf.edu/chimera/download.html
    
    Important: Install Chimera (python 2.7 based) and not ChimeraX (python 3.)
    
    ---
    
    input:  -f          pdb-file with single triple helix
            -dc         contact distance
            -cut-off    cut-off cuts fibril to desired length
           
    output: -cc    crystal-contacts.txt      

    ---

    matrixget:  Calls Chimera via Python 2.7 script from terminal    
    
                -> Gets transformation matrices based on crystal contacts 
                   (e.g. no_cc = 60) to setup skeleton of microfibril
                
    matrixset:  Call Chimera via Python 2.7 script from terminal    
    
                -> set pdb-models based on updated, symmetriyzed transformation
                   matrices.
    swapaa:     Call Chimera via Python 2.7 script from terminal    
    
                -> swaps amino acid defined by user for mutation
    
    """
    def __init__(self,pdb=None):
        self.pdb_file=pdb
        
    def matrixget(self,pdb=None,contact_distance=int,crystalcontacts=str):
        if pdb==None: pdb=self.pdb_file
        path_matrix=os.path.dirname(os.path.realpath(__file__))+'/'
        
        return subprocess.run(
            'chimera --nogui --silent --script "'+str(path_matrix)+'matrixget.py '+
            str(pdb)+'.pdb '+str(contact_distance)+' '+str(crystalcontacts)+'"',
            shell=True,stdout=subprocess.DEVNULL,stderr=subprocess. DEVNULL)        
    
    def matrixset(self,pdb=None,crystalcontacts=str,system_size=int,fibril_length=float):
        if pdb==None: pdb=self.pdb_file
        path_matrix=os.path.dirname(os.path.realpath(__file__))+'/'

        return subprocess.run(
            'chimera --nogui --silent --script "'+str(path_matrix)+'matrixset.py '+
            str(pdb)+'.pdb '+str(crystalcontacts)+' '+str(system_size)+' '+
            str(fibril_length)+'"',shell=True,stdout=subprocess.DEVNULL,stderr=subprocess. DEVNULL)

    def swapaa(self,delete=str,system_type=str):
        path_matrix=os.path.dirname(os.path.realpath(__file__))+'/'

        return subprocess.run(
            'chimera --nogui --silent --script "'+str(path_matrix)+'swapaa.py '+
            str(delete)+' '+str(system_type)+'"',shell=True,stdout=subprocess.DEVNULL,stderr=subprocess. DEVNULL)