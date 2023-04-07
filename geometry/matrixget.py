# Apply Crystal Contacts Command within Chimera To obtain
# Transformation matrices containing information for each TH
import os
import sys
import chimera
from chimera import runCommand as rc
#
path_pdb=str(sys.argv[1])
file_name=str(sys.argv[2])
crystal_dimension=str(sys.argv[3])
crystal_out=str(sys.argv[4])
#
print path_pdb
os.chdir(path_pdb)
#
rc("open "+file_name)
contacts=int(crystal_dimension) 
# Run Crystal Contacts Command
rc("crystalcontacts #0 "+str(contacts)+" copies true schematic false")
# Get Transformation-matrices
rc("matrixget "+crystal_out+".txt")
