import os
import sys
import math
import numpy as np
from chimerax.core.commands import run

pdb_file=str(sys.argv[1])
crystalcontacts_file=str(sys.argv[2])
system_size=int(sys.argv[3])
fibril_length=float(sys.argv[4])

for i in range(0,system_size):
    run(session, "open "+pdb_file+" format pdb")

start_pos=np.array(session.models[0].atoms[0].residue.center)
end_pos=np.array(session.models[0].atoms[-1].residue.center)
center_pos=math.sqrt( (abs(end_pos[2])-abs(start_pos[2]))**2 ) / 2 + start_pos[2]
start_pos[2]=center_pos-.5*fibril_length*10
end_pos[2]=center_pos+.5*fibril_length*10

model_string=" ".join('#'+str(i.id[0]) for i in session.models)
run(session,"open "+crystalcontacts_file+".positions models "+model_string)

for model in range(0,system_size):
    run(session,"save "+str(model)+".pdb #"+str(model+1))
    run(session,"del #"+str(model+1))
    run(session,"open "+str(model)+".pdb")
    os.remove(str(model)+".pdb")

selectAtoms=[]
for model in session.models:
    for residue in model.residues:
        cog=np.mean([atom.coord for atom in residue.atoms],axis=0)
        if cog[2]<=start_pos[2] or cog[2]>=end_pos[2]:
            selectAtoms.append(residue.atomspec)

run(session,"select "+"".join(atom for atom in selectAtoms))
run(session,"del sel")

with open(str(crystalcontacts_file).replace('_opt','')+"_id.txt",'w') as f:
    for model in session.models:
        run(session,"save "+str(model.id[0])+".pdb #"+str(model.id[0]))
        f.write('Model '+str(model.id[0])+'\n')
f.close()
exit()
