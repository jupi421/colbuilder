import subprocess
import os
from colbuilder.geometry import system
from colbuilder.topology import amber, martini                              

def build_martini3(system: system.System,force_field=None,topology_file=None) -> martini.Martini:
    """
    
    build martini 3 force field topology
    
    """
    ff=force_field
    martini_=martini.Martini(system=system,force_field=ff)
    for model_id in system.get_models():
        for connect_id in system.get_model(model_id=model_id).connect:
            martini_.set_pdb(connect_id=connect_id)
            martini_.write_pdb(connect_id=connect_id)

            #subprocess.run(
            #'martinize2 -f order.pdb -sep -merge A,B,C -collagen -ff '+str(ff)+'C'+
            #'-from amber99 -x ')
    
    return martini_
                                 
def build_amber99(system: system.System,force_field=None,topology_file=None) -> amber.Amber:
    """"
    
    builder amber99 force field topology
    
    """
    ff=force_field+'sb-star-ildnp'
    amber_=amber.Amber(system=system,ff=ff)
    try:
        subprocess.run('cp '+str(ff)+'.ff/residuetypes.dat .',shell=True)
        subprocess.run('cp '+str(ff)+'.ff/specbond.dat .',shell=True)
    except:
        print('Error: No force field. Get '+str(ff)+'-force field from colbuilder 1.0.')
        exit()

    print('-- Run pdb2gmx with GROMACS using '+str(ff)+'.ff --')

    for model in system.get_models():
        type_=amber_.merge_pdbs(connect_id=model)

        if type_!=None and os.path.exists(os.getcwd()+'/'+str(type_)+'/'+str(int(model))+'.merge.pdb'): 

            subprocess.run(
            'gmx pdb2gmx -f '+str(type_)+"/"+str(int(model))+'.merge.pdb -ignh '+'-merge all '+
            '-ff '+str(ff)+' -water tip3p -p col_'+str(int(model))+'.top '+
            '-o col_'+str(int(model))+'.gro '+'-i posre_'+str(int(model))+'.itp',shell=True,
            stdout=subprocess.DEVNULL,stderr=subprocess. DEVNULL)

            amber_.write_itp(itp_file='col_'+str(int(model))+'.top')

    return amber_

def build_topology(system: system.System,force_field=None,top_file=None,gro_file=None) -> system.System:
    """
    
    builds topology of a system
    
    """
    if force_field=='amber99':
        ff=force_field+'sb-star-ildnp'
        print('-- Build topology based on '+str(ff)+'.ff --')

        amber_=build_amber99(system=system,force_field=force_field)
        amber_.write_topology(system=system,topology_file=top_file)
        amber_.write_gro(system=system,gro_file=gro_file)
    
    if force_field=='martini3':
        ff=force_field+'C'
        print('-- Build topology based on '+str(ff)+' --')

        martini_=build_martini3(system=system,force_field=force_field)

    return system