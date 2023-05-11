import os
import sys
import math
from chimera import runCommand as rc
from chimera import openModels,selection,Point

pdb_file=str(sys.argv[1])
crystalcontacts_file=str(sys.argv[2])
system_size=int(sys.argv[3])
fibril_length=float(sys.argv[4]) * 10 # from nm to Angstrom



for i in range(0,system_size): 
    rc("open "+pdb_file)

start_pos=openModels.list()[0].atoms[0].coord()
end_pos=openModels.list()[0].atoms[-1].coord()
center_pos=math.sqrt( (abs(end_pos[2])-abs(start_pos[2]))**2 ) / 2 + start_pos[2]
start_pos[2]=center_pos-.5*fibril_length
end_pos[2]=center_pos+.5*fibril_length

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
        if com[2]>=start_pos[2] and com[2]<=end_pos[2]:
            for atom in residue.atoms: 
                selectAtoms.append(atom)

selectFibril=selection.ItemizedSelection()
selectFibril.add(selectAtoms)
selectFibril.addImplied()
selection.mergeCurrent(selection.REPLACE,selectFibril)
rc("sel invert")
rc("del sel")

with open(str(crystalcontacts_file).replace('_opt','')+"_id.txt",'w') as f:
    for model in openModels.list(): 
        rc("write #"+str(model.id)+' '+str(model.id)+".pdb")
        f.write('Model '+str(model.id)+'\n')
f.close()