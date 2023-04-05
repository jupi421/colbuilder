#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:17:32 2023

@author: broszms
"""
import sys
import subprocess
#
def matrixget(pdb,number_crystalcontacts):
    # Call Chimera via Python 2.8 script through Terminal
    # Using a subprocess:
    # -> Matrixget gets the transformation matrices using a user-defined 
    # CrystalContact distance (e.g. 60) as input
    file=str(pdb)
    no_cc=int(number_crystalcontacts)
    subprocess.run('chimera --nogui --silent --script '+'"matrixget.py '+str(file)+' '+str(no_cc)+'"',shell=True)
#
def matrixset(pdb,number_crystalcontacts):
    # Call Chimera via Python 2.8 script through Terminal
    # Using a subprocess:
    # -> Matrixget gets the transformation matrices using a user-defined 
    # CrystalContact distance (e.g. 60) as input
    file=str(pdb)
    no_cc=int(number_crystalcontacts)
    subprocess.run('chimera --nogui --silent --script '+'"matrixget.py '+str(file)+str(no_cc)+'"',shell=True)

if __name__ == '__main__':
    filename='Rat.pdb'#str(sys.argv[1])
    crystalcontacts=2#int(sys.argv[2])
    print('TODO')