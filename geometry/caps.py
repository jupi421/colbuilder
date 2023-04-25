import os
import sys
from pymol import cmd, editor


class Caps:
    """

    Adding Caps to triple helices

    Original version implemented by Dr. Agnieszka Obarska-Kosinska from  
    
    Obarska-Kosinska A, Rennekamp B, Ünal A, Gräter F. 
    ColBuilder: A server to build collagen fibril models.
    Viophys J. 2021 Sep 7;120(17):3544-3549. 
    doi: 10.1016/j.bpj.2021.07.009. 
    Epub 2021 Jul 13. PMID: 34265261; PMCID: PMC8456305.


    """
    def __init__(self,pdb=None):
        self.pdb_file=pdb