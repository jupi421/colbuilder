import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
import pprint as pp
import functions_v4 as func

gro = 'ions.gro'

u = MDAnalysis.Universe(gro)
print (u)

ACE_CA = u.select_atoms('resname ACE and name CH3')
pp.pprint (ACE_CA)

index_file = 'index.ndx'
index_new = 'index_pull.ndx'
data_all, data_array = func.get_data_from_file(index_file)
ctr = 1
ACE_all = '[ ACE_all ] \n'

for i in range(0,len(ACE_CA),3):
    new_line1 = '[ ACE_' +str(ctr) +' ] \n'
    ctr +=1
    new_line2 = str(ACE_CA[i].index + 1) + ' ' + str(ACE_CA[i+1].index +1 ) + ' ' + str(ACE_CA[i+2].index +1) + ' \n'
    ACE_all += str(ACE_CA[i].index + 1) + ' ' + str(ACE_CA[i+1].index +1 ) + ' ' + str(ACE_CA[i+2].index +1) + ' '
    data_all.append(new_line1)
    data_all.append(new_line2)

ACE_all = ACE_all + ' \n'
data_all.append(ACE_all)
    
NME_CA = u.select_atoms('resname NME and name CH3')
NME_all = '[ NME_all ] \n'

ctr = 1
for i in range(0,len(NME_CA),3):
    new_line1 = '[ NME_' +str(ctr) +' ] \n'
    ctr +=1
    new_line2 = str(NME_CA[i].index + 1) + ' ' + str(NME_CA[i+1].index +1 ) + ' ' + str(NME_CA[i+2].index +1) + ' \n'
    NME_all += str(NME_CA[i].index + 1) + ' ' + str(NME_CA[i+1].index +1 ) + ' ' + str(NME_CA[i+2].index +1) + ' '
    data_all.append(new_line1)
    data_all.append(new_line2)

NME_all = NME_all  + ' \n'
data_all.append (NME_all)
    
func.store_linelist_to_file(data_all, index_new)
