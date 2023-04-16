import os
import sys
import chimera
from chimera import runCommand as rc

path_pdb=str(sys.argv[1])
file_name=str(sys.argv[2])
contacts=str(int(str(sys.argv[3])))
crystal_out=str(sys.argv[4])

os.chdir(path_pdb)
rc("open "+file_name)
rc("crystalcontacts #0 "+contacts+" copies true schematic false")
rc("matrixget "+crystal_out+".txt")
