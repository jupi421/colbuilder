# This Script is a first Test Run to Generate the
# Largest Collagen Fibril in the Whole World!
import os
import sys
import chimera
from chimera import runCommand as rc
from chimera import openModels,selection,Point
from subprocess import call
#
#
# Change Directory To File
wd=os.getcwd()
os.chdir(wd)
fileName="Rat.pdb"
#
rc("close session")
rc("open "+fileName)
#
m=openModels.list()
startPos=m[0].atoms[0].coord()
endPos=m[0].atoms[-1].coord()
startPos[2]=startPos[2]-150
endPos[2]=endPos[2]+150
#
# Open Create Crystal Contacts
Contacts=int(60) 
#rc("crystalcontacts #0 "+str(Contacts)+" copies true schematic false")
# Get Trans-Rot-Matrix
#rc("matrixget CrystalContacts.txt")
# Close & Re-Open Session
rc("close session")
# Symmetrize Trans-Rot-Matrix
#call("python crystalSymmetry.py",shell=True)
# Get Length of sym. models
with open('lenModelsSym.txt','r') as f:
    lenModels=[int(l) for l in f][0]
f.close()
#
# Open default model (N=lenModelsMod) times
for i in range(0,lenModels):
    rc("open "+fileName)
# Set Trans-Rot-Matrix with regard to sym.info
rc("matrixset CrystalContactsSym.txt")
#
cnt=1
# Write Trans-Rot models to pdb
for i in openModels.list():
    rc("write #"+str(i.id)+" "+str(i.id+1)+".pdb")
    rc("del #"+str(i.id))
    #
    rc("open "+str(cnt)+".pdb")
    os.remove(str(cnt)+".pdb")
    cnt+=1
#
# select residues within 300 nm distance
selectedAtoms=[]
for i in range(0,lenModels):
    tmp=openModels.list()[i]
    for r in tmp.residues:
        mass=[a.element.mass for a in r.atoms]
        com=Point([a.coord() for a in r.atoms],mass)
        if com[2]>=startPos[2] and com[2]<=endPos[2]:
            for a in r.atoms:
                selectedAtoms.append(a)
#
# Prepare Cutted Fibril: Selection, Invert, Delete
fibrilSelection=selection.ItemizedSelection()
fibrilSelection.add(selectedAtoms)
fibrilSelection.addImplied()
selection.mergeCurrent(selection.REPLACE,fibrilSelection)
rc("sel invert")
rc("del sel")
#
# Save finalized PDBs
for i in openModels.list():
    rc("write #"+str(i.id)+' '+str(i.id+1)+".pdb")
