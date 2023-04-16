import os
import sys
import chimera
from chimera import runCommand as rc
from chimera import openModels,selection,Point

path_pdb=str(sys.argv[1])
file_name=str(sys.argv[2])
crystalcontacts=str(sys.argv[3])
fibril_number=int(sys.argv[4])
cut_off=float(sys.argv[5])

os.chdir(path_pdb)

rc("close")
for i in range(0,fibril_number):
    rc("open "+file_name)

model=openModels.list()
startPos=model[0].atoms[0].coord()
endPos=model[0].atoms[-1].coord()
cut_off=int(cut_off-300)/2
startPos[2]=startPos[2]-cut_off
endPos[2]=endPos[2]+cut_off
rc("matrixset "+crystalcontacts+".txt")

cnt=1
for m in model:
    rc("write #"+str(m.id)+" "+str(m.id+1)+".pdb")
    rc("del #"+str(m.id))
    rc("open "+str(cnt)+".pdb")
    os.remove(str(cnt)+".pdb")
    cnt+=1

selectAtoms=[]
for i in range(0,fibril_number):
    m=model[i]
    for r in m.residues:
        mass=[a.element.mass for a in r.atoms]
        com=Point([a.coord() for a in r.atoms],mass)
        if com[2]>=startPos[2] and com[2]<=endPos[2]:
            for a in r.atoms:
                selectAtoms.append(a)

selectFibril=selection.ItemizedSelection()
selectFibril.add(selectAtoms)
selectFibril.addImplied()
selection.mergeCurrent(selection.REPLACE,selectFibril)
rc("sel invert")
rc("del sel")

for m in model:
    rc("write #"+str(m.id)+' '+str(m.id+1)+".pdb")
