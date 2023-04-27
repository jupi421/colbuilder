"""

Module to build the collagen microfibril from a single collagen triple helix

"""

from colbuilder.geometry import (
    crystal, crystalcontacts, chimera, connect,
    model, system, optimize, caps
)

def build_geometry(path_wd=str,pdb_file=None,contact_distance=float,
                      crystalcontacts_file=str,crystalcontacts_optimize=bool,cut_off=float):
    """
    
    Builds a system of models depending on the input-data from the user.

    i) pdb != None & contact_distance != None & crystalcontacts == None : build system from zero
    ii) pdb != None & contact_distance == None & crystalcontacts != None: build system based on crystalcontacts
    
    """
    if pdb_file==None: print('Error: No pdb-file given to build collagen microfibril')

    print('-- Read crystallographic symmetry from '+str(pdb_file)+'.pdb --')
    crystal_=crystal.Crystal(pdb_file)

    path_pdb_file=path_wd+'/'+pdb_file
    chimera_=chimera.Chimera(path_pdb_file)

    if pdb_file!=None and contact_distance!=None and crystalcontacts_file=='chimera_crystalcontacts':
        """
    
        PART ONE: Generate System of Models based on Contact Distance and PDB-File of Single Molecule
    
        """
        print('-- Build System from '+str(pdb_file)+'.pdb and Contact Distance '+str(contact_distance)+' A --')
        chimera_.matrixget(pdb=path_pdb_file,contact_distance=contact_distance,crystalcontacts=crystalcontacts_file)

        print('-- Write CrystalContacts of Systems to '+str(crystalcontacts_file)+'.txt --')
        crystalcontacts_=crystalcontacts.CrystalContacts(crystalcontacts_file)

        print('-- Build System from '+str(crystalcontacts_file)+'.txt --')
        system_=build_system(crystal=crystal_,crystalcontacts=crystalcontacts_)

        print('-- Connect Models System built --')
        system_,connect_=connect_system(system=system_)

        print('-- Optimize System by Adding Models --')
        system_=optimize.Optimizer(system=system_).run_optimize(system=system_,connect=connect_)

        print('-- Connect Models in Optimized System --')
        system_,connect_=connect_system(system=system_)

        print('-- Write Connected Models from Optimized System --')
        crystalcontacts_.write_contacts(system=system_,crystalcontact_file=crystalcontacts_file+'_connect')

    elif pdb_file!=None and contact_distance==None:
        """
    
        PART TWO: Generate System of Models based user-specific CrystalContacts and PDB-File of Single Molecule
    
        """
        print('-- Build System from PDB structure and user-specified CrystalContacts --')
        crystalcontacts_=crystalcontacts.CrystalContacts(crystalcontacts_file)

        print('-- Build System from '+str(crystalcontacts_file)+'.txt --')
        system_=build_system(crystal=crystal_,crystalcontacts=crystalcontacts_)

        print('-- Connect Models in System built from '+str(crystalcontacts_file)+'.txt --')
        system_,connect_=connect_system(system=system_)

        # TODO: Output connecte models & unconnected ones        
        if crystalcontacts_optimize:
            print('-- Optimize user-specific CrystalContacts by Adding Models --')
            system_=optimize.Optimizer(system=system_).run_optimize(system=system_,connect=connect_)

            print('-- Connect Models in Optimized System built from '+str(crystalcontacts_file)+'.txt --')
            system_,connect_=connect_system(system=system_)

            print('-- Write CrystalContacts Models from Optimized System to '+str(crystalcontacts_file)+'.txt --')
            crystalcontacts_.write_contacts(system=system_,crystalcontact_file=crystalcontacts_file+'_connect')

    else:
        print('Error: Input does not make sense, please provide either Contact Distance or CrystalContacts and not both.')
        return exit()
    
    print('-- Generate Microfibril from '+str(crystalcontacts_file)+' with Chimera --')
    print('-- Please wait, this may take some time ... --')
    chimera_.matrixset(pdb=pdb_file,crystalcontacts=crystalcontacts_file+'_connect',
                       system_size=system_.size_system(system=system_),cut_off=cut_off)   
    
    print('-- Cap N- and C- termini of each Model within Microfibril --')
    caps_=cap_system(system=system_)

    print('-- Merge capped Models to Microfibril --')
    system_=system_.write_system_pdb(system=system_)

    return system_

def build_system(crystal : crystal.Crystal,
                 crystalcontacts : crystalcontacts.CrystalContacts) -> system.System:
    """
    
    Build a system of models
    
    """
    system_=system.System(crystal=crystal,crystalcontacts=crystalcontacts)
    model_t=system_.crystalcontacts.read_t_matrix()
    model_s={ k:system_.crystal.get_s_matrix(t_matrix=model_t[k]) for k in model_t }
    for key_m in model_t:
        model_=model.Model(model_id=key_m,model_t=model_t[key_m],model_s=model_s[key_m])
        system_.add_model(model=model_)
    return system_

def connect_system(system : system.System) -> tuple[system.System, connect.Connect] :
    """
    
    Identifies connections within a build system 
    
    """
    connect_=connect.Connect(system.crystal.pdb_file)
    system_connect=connect_.run_connect(system=system)
    for key_m in system_connect:
        system.get_model(model_id=key_m).add_model_connect(model_connect_id=key_m,model_connect=system_connect[key_m])
    return system,connect_

def cap_system(system : system.System):
    """
    
    Caps each model of system
    
    """
    caps_=caps.Caps(system=system)
    for model_id in range(system.size_system(system=system)):
        caps_.read_residues(pdb_id=int(model_id))
        caps_.add_caps(pdb_id=int(model_id))
    return caps_