
import os
import subprocess

from colbuilder.geometry import system
from topology import amber
                                 
def prepare_pdbs(system=None,force_field=None):
    """"
    
    prepare pdbs for  force field
    
    """
    amber_=amber.Amber(system=system)
    if force_field=='amber99':
        for model in system.get_keys():
            amber_.merge_pdbs(connect_id=model)
                                 
def build_amber99(system=None,force_field=None):
    """"
    
    builder amber99 force field
    
    """
    amber_=amber.Amber(system=system)
    ff=force_field+'sb-star-ildnp'
    try:
        subprocess.run('cp '+str(ff)+'.ff/residuetypes.dat .',shell=True)
        subprocess.run('cp '+str(ff)+'.ff/specbont.dat .',shell=True)
    except:
        print('Error: No force field. Get '+str(ff)+'-force field from colbuilder 1.0.')
        exit()

    for model in system.get_keys():
        amber_.merge_pdbs(connect_id=model)

        subprocess.run(
        'gmx pdb2gmx -f '+str(int(model))+'.merge.pdb -ignh '+'-merge all '+
        '-ff '+str(ff)+' -water tip3p -p col_'+str(int(model))+'.top '+
        '-o col_'+str(int(model))+'.gro '+'-i posre_'+str(int(model))+'.itp',shell=True)

        amber_.
        
def build_topology(system : system.System,force_field=None):
    """
    
    builds topology of a system
    
    """
    prepare_pdbs(system=system,force_field=force_field)
    if force_field=='amber99':
        build_amber99(system=system,force_field=force_field)

