"""

Module to build the collagen microfibril from a single collagen triple helix

"""
import subprocess
from colbuilder.geometry import (
    crystal, crystalcontacts, chimera, crosslink,
    model, system, connect, caps, 
    optimize, mix, mutate,
)

def mutate_geometry(path_wd=None,crystalcontacts_file=None,crystalcontacts_optimize=None,fibril_length=None,
                 pdb_files=None,setup_mutate=None,system=system.System) -> system.System:
    """
    
    A built system is mutated according to user-input parameters.
    Mutation is random, however constrained to ensure at least one crosslink connection
    per side of each triple helix

    correct: --- X-M ---- M-X ----  with M= Mutation, X= Crosslink, --- triple helix
    false:  --- M-M ---- X-X ---- 

    """



    return print('mutate')


def mix_geometry(path_wd=None,crystalcontacts_file=None,crystalcontacts_optimize=None,fibril_length=None,
                 pdb_files=None,setup_mix=None,system=system.System,pdb_out=str) -> system.System:
    """
    
    A built system is mixed according to user-input parameters.
    Mixing is random, however constrained due to ratio of mixed models.
    
    """
    mix_setup={ idx.split(':')[0]:idx.split(':')[1] for idx in setup_mix }
    mix_pdb=dict(zip(mix_setup.keys(),pdb_files))

    print('-- Prepare Mix '+str(setup_mix)+' --')
    if crystalcontacts_optimize: crystalcontacts_file=crystalcontacts_file+'_opt'

    mix_setup={ idx.split(':')[0]:idx.split(':')[1] for idx in setup_mix }
    mix_pdb=dict(zip(mix_setup.keys(),pdb_files))

    chimera_=chimera.Chimera(path_wd+'/'+pdb_files[0])
    for key in list(mix_pdb.keys())[1:]:
        subprocess.run('mkdir '+path_wd+'/'+key,shell=True)

        print('-- Generate System for '+str(key)+' with '+str(mix_pdb[key])+' --')
        chimera_.matrixset(pdb=mix_pdb[key],crystalcontacts=crystalcontacts_file,
                    system_size=system.get_size(system=system),fibril_length=fibril_length)
        
        print('-- Add Caps --')
        cap_system(system=system)

        subprocess.run('mv *.caps.pdb '+key,shell=True)

    system_=mix.Mix(setup=mix_setup,system=system).add_mix()

    print('-- Mix System --')
    connect_file=crystalcontacts_file.replace('opt','')+'_connect_mix'
    connect.Connect(pdb=pdb_files[0]).write_mix_connect(system=system_,connect_file=connect_file)

    system_.write_pdb(pdb_out=pdb_out)

    return system_

def build_geometry(path_wd=str,pdb_file=None,contact_distance=float,crystalcontacts_file=str,
                   crystalcontacts_optimize=bool,fibril_length=float,pdb_out=str)  -> system.System:
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

    print('-- Generate System from '+str(crystalcontacts_file)+'_opt --')
    print('-- Please wait, this may take some time ... --')
    
    chimera_.matrixset(pdb=pdb_file,crystalcontacts=crystalcontacts_file+'_opt',
                       system_size=system_.get_size(system=system_),fibril_length=fibril_length)   
    
    print('-- Update System after Matrixset --')
    system_=matrixset_system(system_=system_,crystalcontacts_=crystalcontacts_,
                             crystalcontacts_file=crystalcontacts_file+'_opt')

    print('-- Add Caps --')
    subprocess.run('mkdir '+path_wd+'/'+system_.get_model(model_id=0.0).type,shell=True)
    cap_system(system=system_)
    subprocess.run('mv *.caps.pdb '+system_.get_model(model_id=0.0).type,shell=True)

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
        model_=model.Model(id=key_m,transformation=transformation[key_m],unit_cell=unit_cell[key_m],pdb_file=crystal.pdb_file)
        system_.add_model(model=model_)
    return system_

def connect_system(system : system.System) -> tuple[system.System, connect.Connect]:
    """
    
    Identifies connections within a build system 
    
    """
    connect_=connect.Connect(system=system)
    system_connect=connect_.run_connect(system=system)
    for key_m in system_connect:
        system.get_model(model_id=key_m).add_connect(connect_id=key_m,connect=system_connect[key_m])
    return system,connect_

def cap_system(system : system.System) -> caps.Caps:
    """
    
    Cap each model of system
    
    """
    caps_=caps.Caps(system=system)
    for idx in system.get_keys():
        caps_.read_residues(pdb_id=int(idx))
        caps_.add_caps(pdb_id=int(idx))
    return caps_

def matrixset_system(system_ : system.System,crystalcontacts_ : crystalcontacts.CrystalContacts,
                     crystalcontacts_file=str) -> system.System:
    """
    
    set system after the matrixset command in chimera
    
    """
    contacts=crystalcontacts_.get_models(crystalcontact_file=crystalcontacts_file)
    for model_id in system_.get_keys():
        if model_id not in contacts: system_.delete_model(model_id=model_id)
    return system_

def write_system(system_ : system.System,crystalcontacts_ : crystalcontacts.CrystalContacts,
                 connect_ : connect.Connect):
    """
    
    Wrtie crystalcontacts and model connection within system to txt-file
    
    """
    print('-- Write CrystalContacts '+str(crystalcontacts_.crystalcontact_file)+'_opt --')
    crystalcontacts_.write_crystalcontacts(system=system_,crystalcontact_file=crystalcontacts_.crystalcontact_file+'_opt')
    
    print('-- Write Connected Models '+str(crystalcontacts_.crystalcontact_file)+'_connect --')
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
    chimera_.matrixget(pdb=path_pdb_file,contact_distance=contact_distance,
                       crystalcontacts=crystalcontacts_file)

    print('-- Write CrystalContacts to '+str(crystalcontacts_file)+' --')
    crystalcontacts_=crystalcontacts.CrystalContacts(crystalcontacts_file)

    print('-- Build System --')
    system_=build_system(crystal=crystal_,crystalcontacts=crystalcontacts_)
        
    print('-- Connect System --')
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
    
    print('-- Build System --')
    system_=build_system(crystal=crystal_,crystalcontacts=crystalcontacts_)

    print('-- Connect System --')
    system_,connect_=connect_system(system=system_)
   
    if crystalcontacts_optimize:
        print('-- Optimize System --')
        system_=optimize.Optimizer(system=system_).run_optimize(system=system_,connect=connect_)
        system_,connect_=connect_system(system=system_)

    return system_,crystalcontacts_,connect_