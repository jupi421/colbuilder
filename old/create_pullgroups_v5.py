import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
import pprint as pp
import functions_v4 as func


###create pull groups in indexfile


topfile = 'topol.top'

top_all, top_array = func.get_data_from_file(topfile)

index_file = 'index.ndx'
index_new2 = 'index_pull.ndx'

data_all, data_array = func.get_data_from_file(index_file)

ACE_all = '[ ACE_all ] \n'
NME_all = '[ NME_all ] \n'

residues = []
res_ctr = 0
res_ACE_NME = []
triple_helix_ctr1 = 0
triple_helix_ctr2 = 0
ACE_ctr = 1
ACE_term_ctr = 1
NME_ctr = 1
NME_term_ctr = 1
ACE_group = ''
NME_group = ''

for j in range(len(top_all)):
    if len(top_array[j]) < 4:
        continue
    if 'residue' in top_array[j][1]:
        residue = (top_array[j][2], top_array[j][3])
        residues.append(residue)
    if (('ACE' in top_array[j][3]) and ('CH3' in top_array[j][4])):     #take CH3-atom          
        ACE_all += str(top_array[j][0]) + ' ' 
        triple_helix_ctr1 += 1 
        ACE_group += str(top_array[j][0]) + ' '
        residue = (top_array[j][2], top_array[j][3])
        res_ctr += 1
        res_ACE_NME.append(residue)
        if triple_helix_ctr1 == 3:
            if res_ACE_NME[res_ctr-2][1] == 'NME':
                new_line = '[ ACE_' +str(ACE_ctr) +' ] \n' + ACE_group     + '\n'         
                data_all.append(new_line)
                ACE_ctr += 1
                triple_helix_ctr1 = 0
                ACE_group = ''
            else:
                new_line = '[ ACE_term_' +str(ACE_term_ctr) +' ] \n' + ACE_group     + '\n'         
                data_all.append(new_line)
                ACE_term_ctr += 1
                triple_helix_ctr1 = 0
                ACE_group = ''

    if 'NME' in top_array[j][3] and 'CH3' in top_array[j][4]:               
        NME_all += str(top_array[j][0]) + ' '   
        triple_helix_ctr2 += 1 
        NME_group += str(top_array[j][0]) + ' '
        residue = (top_array[j][2], top_array[j][3])
        res_ACE_NME.append(residue)
        res_ctr += 1
        if triple_helix_ctr2 == 3:
            if res_ACE_NME[res_ctr-2][1] == 'ACE':
                new_line = '[ NME_' +str(NME_ctr) +' ] \n' + NME_group     + '\n'  
                data_all.append(new_line)
                NME_ctr += 1
                triple_helix_ctr2 = 0
                NME_group = ''
            else:
                new_line = '[ NME_term_' +str(NME_term_ctr) +' ] \n' + NME_group     + '\n'         
                data_all.append(new_line)
                NME_term_ctr += 1
                triple_helix_ctr2 = 0
                NME_group = ''           


            
#pp.pprint (res_ACE_NME)
ACE_all += '\n'
NME_all += '\n'
data_all.append(ACE_all)
data_all.append(NME_all) 
func.store_linelist_to_file(data_all, index_new2)
    


