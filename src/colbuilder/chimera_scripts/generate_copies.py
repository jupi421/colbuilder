import chimera
from chimera import runCommand
import os

input_pdb = os.environ.get('INPUT_PDB')
model = chimera.openModels.open(input_pdb)[0]

runCommand("crystalcontacts #0 1.5 copies true schematic false")

models = chimera.openModels.list(modelTypes=[chimera.Molecule])
generated_pdbs = []

for i, m in enumerate(models):
    filename = "copy_{}.pdb".format(i)
    runCommand("write format pdb #{} {}".format(m.id, filename))
    generated_pdbs.append(filename)

with open("generated_pdbs.txt", "w") as f:
    for pdb in generated_pdbs:
        f.write(pdb + "\n")

print("Generated {} PDB files.".format(len(generated_pdbs)))

runCommand("close #")