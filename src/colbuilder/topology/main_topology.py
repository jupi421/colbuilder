"""

module to build the topology for the system

"""
import subprocess
from colbuilder.geometry import ( system,
)
from colbuilder.topology import ( merge
)
                                 
def prepare_pdbs(system=None,force_field=None):
    """"
    
    prepare pdbs for  force field
    
    """
    merge_=merge.Merge(system=system)
    if force_field=='amber99':
        for model in system.get_keys():
            merge_.merge_pdbs(connect_id=model)
                                 
def build_amber99(system=None):
    """"
    
    builder amber99 force field
    
    """
    merge_=merge.Merge(system=system)
    for model in system.get_keys():
        merge_.merge_pdbs(connect_id=model)
        subprocess.run(
        'printf '"1\n6\n"'gmx pdb2gmx -f '+str(int(model))+'.merge.pdb -ignh'+ 
        '-merge all'+'-p col_'+str(model)+'.top'+'-o col_'+str(model)+'.gro'+
        '-p posre_'+str(model)+'.itp',shell=True)

def build_topology(system : system.System,force_field=None):
    """
    
    builds topology of a system
    
    """
    prepare_pdbs(system=system,force_field=force_field)
    if force_field=='amber99':
        build_amber99(system=system)

