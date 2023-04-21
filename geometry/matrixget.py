import os
import sys
import chimera
from chimera import runCommand as rc


file_name=str(sys.argv[1])
d_contact=str(int(str(sys.argv[2])))
crystal_out=str(sys.argv[3])

rc("open "+file_name)
rc("crystalcontacts #0 "+d_contact+" copies true schematic false")
rc("matrixget "+crystal_out+".txt")
