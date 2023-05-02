"""

Module to generate topology for the collagen microfibril 

"""
from colbuilder.geometry import ( 
    system, connect
) 

from colbuilder.topology import ( 
    topology
) 

def build_cg_topology(system_ : system.System,pdb_file=None,connect=None):
    """
    
    Function builds coarse-grained topology from microfibril
    
    """
    if connect_==None: connect_=connect.Connect(pdb_file)
    if system_==None: system_=system.System()
    print(system_)
    if system_==None:
        system_.read_pdb_system(system_pdb_file=pdb_file)
        martini_=topology.Martini(system_pdb_file=pdb_file)
    print('to be done')

    return 