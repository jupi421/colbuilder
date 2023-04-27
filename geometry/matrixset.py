import os
import sys
from chimera import runCommand as rc
from chimera import openModels,selection,Point

pdb_file=str(sys.argv[1])
crystalcontacts_file=str(sys.argv[2])
system_size=int(sys.argv[3])
cut_off=float(sys.argv[4])

for i in range(0,system_size): 
    rc("open "+pdb_file)

startPos=openModels.list()[0].atoms[0].coord()
endPos=openModels.list()[0].atoms[-1].coord()
cut_off=int(cut_off-300)/2
startPos[2]=startPos[2]-cut_off
endPos[2]=endPos[2]+cut_off

rc("matrixset "+crystalcontacts_file+".txt")

cnt=0
for model in openModels.list():
    rc("write #"+str(model.id)+" "+str(model.id)+".pdb")
    rc("del #"+str(model.id))
    rc("open "+str(cnt)+".pdb")
    os.remove(str(cnt)+".pdb")
    cnt+=1

selectAtoms=[]
for model in openModels.list():
    for residue in model.residues:
        mass=[atom.element.mass for atom in residue.atoms]
        com=Point([atom.coord() for atom in residue.atoms],mass)
        if com[2]>=startPos[2] and com[2]<=endPos[2]:
            for atom in residue.atoms: 
                selectAtoms.append(atom)

selectFibril=selection.ItemizedSelection()
selectFibril.add(selectAtoms)
selectFibril.addImplied()
selection.mergeCurrent(selection.REPLACE,selectFibril)
rc("sel invert")
rc("del sel")

for model in openModels.list(): 
    rc("write #"+str(model.id)+' '+str(model.id)+".pdb")