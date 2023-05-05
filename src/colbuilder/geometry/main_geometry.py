"""

Module to build the collagen microfibril from a single collagen triple helix

"""
import subprocess
from colbuilder.geometry import (
    crystal, crystalcontacts, chimera, connect,
    model, system, optimize, caps, mix, mut
)

def mutate_geometry(path_wd=None,crystalcontacts_file=None,crystalcontacts_optimize=None,cut_off=None,
                 pdb_files=None,setup_mutate=None,system=system.System) -> system.System:
    """
    
    A built system is mutated according to user-input parameters.
    Mutation is random, however constrained to ensure at least one crosslink connection
    per side of each triple helix

    correct: --- X-M ---- M-X ----  with M= Mutation, X= Crosslink, --- triple helix
    false:  --- M-M ---- X-X ---- 

    """



    return print('mutate')


def mix_geometry(path_wd=None,crystalcontacts_file=None,crystalcontacts_optimize=None,cut_off=None,
                 pdb_files=None,setup_mix=None,system=system.System,pdb_out=str) -> system.System:
    """
    
    A built system is mixed according to user-input parameters.
    Mixing is random, however constrained due to ratio of mixed models.
    
    """
    mix_setup={ idx.split(':')[0]:idx.split(':')[1] for idx in setup_mix }
    mix_pdb=dict(zip(mix_setup.keys(),pdb_files))

    # Process built System
    subprocess.run('mkdir '+path_wd+'/'+str(list(mix_setup.keys())[0]),shell=True)
    subprocess.run('mv *.caps.pdb '+str(list(mix_setup.keys())[0]),shell=True)

    print('-- Prepare Mix '+str(setup_mix)+' --')
    if crystalcontacts_optimize: crystalcontacts_file=crystalcontacts_file+'_opt'

    mix_setup={ idx.split(':')[0]:idx.split(':')[1] for idx in setup_mix }
    mix_pdb=dict(zip(mix_setup.keys(),pdb_files))

    chimera_=chimera.Chimera(path_wd+'/'+pdb_files[0])
    for key in list(mix_pdb.keys())[1:]:
        subprocess.run('mkdir '+path_wd+'/'+key,shell=True)

        print('-- Generate Structure for '+str(key)+' with '+str(mix_pdb[key])+' --')
        chimera_.matrixset(pdb=mix_pdb[key],crystalcontacts=crystalcontacts_file,
                    system_size=system.get_size(system=system),cut_off=cut_off)
        
        print('-- Add Caps --')
        cap_system(system=system)

        subprocess.run('mv *.caps.pdb '+key,shell=True)

    system_=mix.Mix(setup=mix_setup,system=system).add_mix()

    print('-- Write Mix Connect --')
    connect_file=crystalcontacts_file.replace('opt','')+'_connect_mix'
    connect.Connect(pdb=pdb_files[0]).write_mix_connect(system=system_,connect_file=connect_file)

    print('-- Write Mix System --')
    system_.write_mix_pdb(pdb_out=pdb_out)

    return system_

def build_geometry(path_wd=str,pdb_file=None,contact_distance=float,crystalcontacts_file=str,
                   crystalcontacts_optimize=bool,cut_off=float,pdb_out=str)  -> system.System:
    """
    
    Builds a system of models from input parameters
    
    """
    if pdb_file==None: print('Error: No pdb-file given to build collagen microfibril')

    print('-- Read crystallographic symmetry from '+str(pdb_file)+'.pdb --')
    crystal_=crystal.Crystal(pdb_file)

    path_pdb_file=path_wd+'/'+pdb_file
    chimera_=chimera.Chimera(path_pdb_file)

    if pdb_file!=None and contact_distance!=None and crystalcontacts_file=='crystalcontacts':
        """

        Generate System of Models from Contact Distance and PDB-File of Single Molecule

        """
        system_,crystalcontacts_,connect_=build_from_contactdistance(path_wd=path_wd,pdb_file=pdb_file,contact_distance=contact_distance,
                                            crystalcontacts_file=crystalcontacts_file,chimera_=chimera_,crystal_=crystal_)

    elif pdb_file!=None and contact_distance==None:
        """

        Generate System of Models from CrystalContacts and PDB-File of Single Molecule

        """
        system_,crystalcontacts_,connect_=build_from_crystalcontacts(crystalcontacts_file=crystalcontacts_file,
                                            crystal_=crystal_,crystalcontacts_optimize=crystalcontacts_optimize)
    else:
        print('Error: Input does not make sense, please provide either Contact Distance or CrystalContacts and not both.')
        return exit()
    
    write_system(system_=system_,crystalcontacts_=crystalcontacts_,connect_=connect_)

    print('-- Generate System from '+str(crystalcontacts_file)+' --')
    print('-- Please wait, this may take some time ... --')
    chimera_.matrixset(pdb=pdb_file,crystalcontacts=crystalcontacts_file+'_opt',
                       system_size=system_.get_size(system=system_),cut_off=cut_off)   
    
    print('-- Add Caps --')
    cap_system(system=system_)

    print('-- Merge Models to Microfibril --')
    system_.write_pdb(pdb_out=pdb_out)

    return system_

def build_system(crystal : crystal.Crystal,
                 crystalcontacts : crystalcontacts.CrystalContacts) -> system.System:
    """
    
    Build a system of models
    
    """
    system_=system.System(crystal=crystal,crystalcontacts=crystalcontacts)
    transformation=system_.crystalcontacts.read_t_matrix()
    unit_cell={ k:system_.crystal.get_s_matrix(t_matrix=transformation[k]) for k in transformation }
    for key_m in transformation:
        model_=model.Model(id=key_m,transformation=transformation[key_m],unit_cell=unit_cell[key_m])
        system_.add_model(model=model_)
    return system_

def connect_system(system : system.System) -> tuple[system.System, connect.Connect]:
    """
    
    Identifies connections within a build system 
    
    """
    connect_=connect.Connect(system.crystal.pdb_file)
    system_connect=connect_.run_connect(system=system)
    for key_m in system_connect:
        system.get_model(model_id=key_m).add_connect(connect_id=key_m,connect=system_connect[key_m])
    return system,connect_

def cap_system(system : system.System) -> caps.Caps:
    """
    
    Cap each model of system
    
    """
    caps_=caps.Caps(system=system)
    for model_id in range(system.size):
        caps_.read_residues(pdb_id=int(model_id))
        caps_.add_caps(pdb_id=int(model_id))
    return caps_

def write_system(system_ : system.System,crystalcontacts_ : crystalcontacts.CrystalContacts,connect_ : connect.Connect):
    """
    
    Wrtie crystalcontacts and model connection within system to txt-file
    
    """
    print('-- Write CrystalContacts '+str(crystalcontacts_.crystalcontact_file)+'_opt.txt --')
    crystalcontacts_.write_crystalcontacts(system=system_,crystalcontact_file=crystalcontacts_.crystalcontact_file+'_opt')
    
    print('-- Write Connected Models '+str(crystalcontacts_.crystalcontact_file)+'_connect.txt --')
    connect_.write_connect(system=system_,connect_file=crystalcontacts_.crystalcontact_file+'_connect')

def build_from_contactdistance(path_wd=str,pdb_file=None,contact_distance=float,
                               crystalcontacts_file=str,chimera_=chimera.Chimera,
                               crystal_=crystal.Crystal) -> tuple[system.System, crystalcontacts.CrystalContacts, connect.Connect]:
    """

    Generate System of Models based on Contact Distance and PDB-File of Single Molecule
    
    """
    path_pdb_file=path_wd+'/'+pdb_file

    print('-- Build System Contact Distance '+str(contact_distance)+' A --')
    print('-- Please wait, this may take some time ... --')
    chimera_.matrixget(pdb=path_pdb_file,contact_distance=contact_distance,crystalcontacts=crystalcontacts_file)

    print('-- Write CrystalContacts '+str(crystalcontacts_file)+'.txt --')
    crystalcontacts_=crystalcontacts.CrystalContacts(crystalcontacts_file)

    print('-- Build System --')
    system_=build_system(crystal=crystal_,crystalcontacts=crystalcontacts_)
        
    print('-- Connect Models --')
    system_,connect_=connect_system(system=system_)

    print('-- Optimize System --')
    system_=optimize.Optimizer(system=system_).run_optimize(system=system_,connect=connect_)
    system_,connect_=connect_system(system=system_)

    return system_,crystalcontacts_,connect_

def build_from_crystalcontacts(crystalcontacts_file=str,crystal_=crystal.Crystal,
                               crystalcontacts_optimize=bool) -> tuple[system.System, crystalcontacts.CrystalContacts, connect.Connect]:
    """
    
    Generate System of Models based user-specific CrystalContacts and PDB-File of Single Molecule
    
    """
    crystalcontacts_=crystalcontacts.CrystalContacts(crystalcontacts_file)
    
    print('-- Build System from CrystalContacts --')
    system_=build_system(crystal=crystal_,crystalcontacts=crystalcontacts_)

    print('-- Connect Models --')
    system_,connect_=connect_system(system=system_)
   
    if crystalcontacts_optimize:
        print('-- Optimize System --')
        system_=optimize.Optimizer(system=system_).run_optimize(system=system_,connect=connect_)
        system_,connect_=connect_system(system=system_)

    return system_,crystalcontacts_,connect_