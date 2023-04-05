# This Script is a first Test Run to Generate the
# Largest Collagen Fibril in the Whole World!
import os
import sys
import chimera
from chimera import runCommand as rc
#
wd=os.getcwd()
os.chdir(wd)
#
fileName=str(sys.argv[1])
crystalDimension=str(sys.argv[2])
#
rc("open "+fileName)
Contacts=int(crystalDimension) 
# Run Crystal Contacts Command
rc("crystalcontacts #0 "+str(Contacts)+" copies true schematic false")
# Get Transformation-matrices
rc("matrixget CrystalContacts.txt")
