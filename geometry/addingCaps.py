import os
import sys
from pymol import cmd, editor


def is_coord_line(line):
    return line[0:6] in ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')

def fix_TER(lines):
    prev = None
    out = []
    for line in lines:
        if line[:3] == 'TER':
            line = prev[:26] + '\n'
            line = line[:12] + '    ' + line[16:]
            line = 'TER   ' + line[6:]
        else:
            prev = line
        out.append(line)

    return out




chainA_resi = []
chainB_resi = []
chainC_resi = []
N_ter_resi = []
C_ter_resi = []


pdb_filename = sys.argv[1]

with open(pdb_filename) as f:
    for line in f:
        if is_coord_line(line):
            
            cur_resi_nr = int(line[22:26])
            cur_resi_chain_id = line[21]
            cur_resi_atom_name = line[13:15]
            
            if cur_resi_atom_name == 'CA' and cur_resi_chain_id == 'A':
                chainA_resi.append(int(line[22:26]))
            
            if cur_resi_atom_name == 'CA' and cur_resi_chain_id == 'B':
                chainB_resi.append(int(line[22:26]))

            if cur_resi_atom_name == 'CA' and cur_resi_chain_id == 'C':
                chainC_resi.append(int(line[22:26]))
        
                
edit_lines_N = []
edit_lines_N.append('resi ' + str(chainA_resi[0]) +' and chain A and name N')
edit_lines_N.append('resi ' + str(chainB_resi[0]) +' and chain B and name N')
edit_lines_N.append('resi ' + str(chainC_resi[0]) +' and chain C and name N')

edit_lines_C = []
edit_lines_C.append('resi ' + str(chainA_resi[-1]) + ' and chain A and name C')
edit_lines_C.append('resi ' + str(chainB_resi[-1]) + ' and chain B and name C')
edit_lines_C.append('resi ' + str(chainC_resi[-1]) + ' and chain C and name C')

cmd.load(pdb_filename)
#cmd.remove('name OXT')

for line_N in edit_lines_N:
    if int(line_N.split(' ')[1])==1:
        continue
    else:
        print(line_N)
        cmd.edit(line_N)
        editor.attach_amino_acid("pk1", 'ace',ss=0)

for line_C in edit_lines_C:
    if int(line_C[5:9])==1056 and line_C[20]!='B':
        continue
    elif line_C[20]=='B' and int(line_C[5:9])==1040:
        continue
    else:
        print(line_C)
        cmd.edit(line_C)
        editor.attach_amino_acid("pk1", 'nme',ss=0)

tmp_filename = pdb_filename.replace('.pdb', '.tmp.pdb') 
cmd.save(tmp_filename)

tmp_lines = []

with open(tmp_filename) as f:
    for line in f:
        
        if line[0:3] != 'TER':
            tmp_lines.append(line)

out_lines = []

prev_resi = None
prev_line = ''
for line in tmp_lines:
    if is_coord_line(line):
        cur_resi = int(line[22:26])
        if (prev_resi is not None) and (cur_resi - prev_resi not in (1,0)) and not prev_line.startswith('TER'):
            out_lines.append('TER\n')
        prev_resi = cur_resi
    prev_line = line
    out_lines.append(line)

out_filename = pdb_filename.replace('.pdb', '.caps.pdb') 

with open(out_filename,'w') as f:
    f.write(''.join(out_lines))

cmd.quit
