import os
import sys
import math
from chimera import runCommand as rc 
from chimera import openModels

pdb_file=str(sys.argv[1])
d_contact=str(int(str(sys.argv[2])))
crystalcontacts_file=str(sys.argv[3])

rc("open "+pdb_file)
start_pos=openModels.list()[0].atoms[0].coord()
end_pos=openModels.list()[0].atoms[-1].coord()
center_pos=math.sqrt( (abs(end_pos[2])-abs(start_pos[2]))**2 ) / 2 + start_pos[2]

if center_pos>0:
    shift_pos=str(int(-center_pos))
else:
    shift_pos=str(int(center_pos))

rc("move 0,0,"+shift_pos+" mod #0 ")
rc("write #0 "+'tmp_move.pdb')
rc("del #0 ")
rc("open tmp_move.pdb")
os.remove('tmp_move.pdb')
rc("crystalcontacts #0 "+d_contact+" copies true schematic false")
rc("matrixget "+crystalcontacts_file+".txt")