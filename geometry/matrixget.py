import sys
from chimera import runCommand as rc

pdb_file=str(sys.argv[1])
d_contact=str(int(str(sys.argv[2])))
crystalcontacts_file=str(sys.argv[3])

rc("open "+pdb_file)
rc("crystalcontacts #0 "+d_contact+" copies true schematic false")
rc("matrixget "+crystalcontacts_file+".txt")