import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
import pprint as pp
import functions_v4 as func


###create pull groups in indexfile


pdbfile = 'RatCross.pdb'

pdb_all, pdb_array = func.get_data_from_file(pdbfile)

index_file = 'index.ndx'
index_new2 = 'index_pull2.ndx'

data_all, data_array = func.get_data_from_file(index_file)

ACE_all = '[ ACE_all ] \n'
NME_all = '[ NME_all ] \n'

residues = []
res_ctr = 0
res_ACE_NME = []
triple_helix_ctr1 = 0
triple_helix_ctr2 = 0
ACE_ctr = 1
NME_ctr = 1
ACE_group = ''
NME_group = ''

for j in range(len(pdb_all)):
    if len(pdb_array[j]) < 4:
        continue
    if (('ACE' in pdb_array[j][3]) and ('CH3' in pdb_array[j][2])):     #take CH3-atom          
        ACE_all += str(pdb_array[j][1]) + ' ' 
        triple_helix_ctr1 += 1 
        ACE_group += str(pdb_array[j][1]) + ' '

        res_ctr += 1
        if triple_helix_ctr1 == 3:
            new_line = '[ ACE_' +str(ACE_ctr) +' ] \n' + ACE_group     + '\n'         
            data_all.append(new_line)
            ACE_ctr += 1
            triple_helix_ctr1 = 0
            ACE_group = ''
    if 'NME' in pdb_array[j][3] and 'CH3' in pdb_array[j][2]:               
        NME_all += str(pdb_array[j][1]) + ' '   
        triple_helix_ctr2 += 1 
        NME_group += str(pdb_array[j][1]) + ' '
        res_ctr += 1
        if triple_helix_ctr2 == 3:
            new_line = '[ NME_' +str(NME_ctr) +' ] \n' + NME_group     + '\n'  
            data_all.append(new_line)
            NME_ctr += 1
            triple_helix_ctr2 = 0
            NME_group = ''
            
#pp.pprint (res_ACE_NME)
ACE_all += '\n'
NME_all += '\n'
data_all.append(ACE_all)
data_all.append(NME_all) 
func.store_linelist_to_file(data_all, index_new2)
    


