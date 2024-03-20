import subprocess
from colbuilder.geometry import (
    crystal, crystalcontacts, chimera, model, system, connect, caps, 
    optimize, mix, fibril, replace
)

def build_geometry(path_wd=str,pdb_file=None,contact_distance=float,crystalcontacts_file=str,connect_file=str,
                   crystalcontacts_optimize=bool,solution_space=[],fibril_length=float,geometry=str,pdb_out=str) -> system.System:
    """
    
    build system of models from input
    
    """
    if pdb_file==None: print('Error: No pdb-file given to build collagen microfibril')

    print('-- Read crystallographic symmetry from '+str(pdb_file)+'.pdb --')
    crystal_=crystal.Crystal(pdb_file)
    crystal_.translate_crystal(pdb=pdb_file,translate=[0,0,4000])

    path_pdb_file=path_wd+'/'+pdb_file
    chimera_=chimera.Chimera(path_pdb_file)

    if pdb_file!=None and contact_distance!=None and crystalcontacts_file==None:
        system_,crystalcontacts_,connect_=build_from_contactdistance(path_wd=path_wd,pdb_file=pdb_file,
                                            contact_distance=contact_distance,solution_space=solution_space,
                                            crystalcontacts_file=crystalcontacts_file,chimera=chimera_,crystal=crystal_)

    elif pdb_file!=None and contact_distance==None and crystalcontacts_file!=None:
        system_,crystalcontacts_,connect_=build_from_crystalcontacts(crystalcontacts_file=crystalcontacts_file,
                                            solution_space=solution_space,crystal=crystal_,connect_file=connect_file,
                                            crystalcontacts_optimize=crystalcontacts_optimize)
    
    elif pdb_file!=None and contact_distance==0 and crystalcontacts_file==None:
        system_,crystalcontacts_,connect_=build_from_pdb(path_wd=path_wd,pdb_file=pdb_file,contact_distance=contact_distance,
                               crystalcontacts_file=crystalcontacts_file,chimera=chimera_,crystal=crystal_)    
    
    elif geometry==True:
        print('Error: Please provide Contact Distance or CrystalContacts to generate the microfibril OR \n'+
              'contact distance = 0 to generate the triple helix topology. NOTE: Your contacts file should not be named \n'+
              'crystalcontacts_from_colbuilder.txt, since this is the default one. Please rename, such that message omits.')
        
        return exit()
    
    if geometry==False: 
        print('-- Set -geometry flag to generate microfibrillar structure pdb file --')
        return system_

    print('-- Write '+str(crystalcontacts_.crystalcontacts_file)+' --')       
    crystalcontacts_.write_crystalcontacts(system=system_,crystalcontacts_file=crystalcontacts_.crystalcontacts_file)

    print('-- Generate system from '+str(crystalcontacts_.crystalcontacts_file)+' --')
    print('-- Please wait, this may take some time ... --')

    chimera_.matrixset(pdb=pdb_file,crystalcontacts=crystalcontacts_.crystalcontacts_file,
                    system_size=system_.get_size(system=system_),fibril_length=fibril_length)   
    
    print('-- Cut system to '+str(fibril_length)+' nm --')
    system_=matrixset_system(system=system_,crystalcontacts_file=crystalcontacts_.crystalcontacts_file)

    print('-- Write '+str(connect_.connect_file)+' --')
    connect_.write_connect(system=system_,connect_file=connect_.connect_file)

    print('-- Add caps --')
    subprocess.run('rm -r '+path_wd+'/'+system_.get_model(model_id=0.0).type,shell=True)
    subprocess.run('mkdir '+path_wd+'/'+system_.get_model(model_id=0.0).type,shell=True)
    cap_system(system=system_,crosslink_type=system_.get_model(model_id=0.0).type)

    print('-- Write '+str(pdb_out)+' --')
    system_.write_pdb(pdb_out=pdb_out,fibril_length=fibril_length)

    return system_

def replace_geometry(path_wd=str,ratio_replace=None,system=system.System,
                    fibril_length=float,pdb_out=str,replace_file=False) -> system.System:
    """
    
    replace crosslinks within microfibril to reduce the overall number of crosslinks

    """
    if replace_file==False:
        print('-- Replace '+str(ratio_replace)+'%'+' of crosslinks from microfibril --')
        replace_=replace.Replace(ratio_replace=ratio_replace,system=system,fibril_length=fibril_length)
        system_=replace_.run_replace(system=system,ratio_replace=ratio_replace)

        replace_.write_replace(system=system_,file='replace')  
        replace_file='replace'

    else:
        system_=system

    print('-- Please wait, this may take some time ... --')
    chimera_=chimera.Chimera(path_wd+'/'+system_.crystal.pdb_file)
    chimera_.swapaa(replace=replace_file,system_type=system_.get_model(model_id=0.0).type)

    system_.write_pdb(pdb_out=pdb_out,fibril_length=fibril_length)

    return system_

def mix_geometry(path_wd=str,fibril_length=float,connect_file=None,
                 pdb_files=[],ratio_mix=None,system=system.System,pdb_out=str) -> system.System:
    """
    
    mix differnetly crosslinked models/ triple helices within the system
    according to user-input parameters:
    1) Provide a connect-file with certain crosslink specifications
    (e.g., 1.caps 3.caps ; T) together with the triple helices files (e.g., Rat-T.pdb Rat-D.pdb)
    2) Provide ratio of mixtures between crosslinked triple helices (e.g., T:70 D:30) and 
    the respective files (see 1))
    
    """
    print('-- Prepare mix setup --')   
    if ratio_mix!=None and connect_file==None:
        
        mix_setup={ idx.split(':')[0]:idx.split(':')[1] for idx in ratio_mix }
        mix_pdb=dict(zip(mix_setup.keys(),pdb_files))

        
        system_size=system.get_size(system=system)
        connect_file='connect_from_colbuilder'

        for key in list(mix_pdb.keys())[1:]:
            
            crystal.Crystal(pdb=mix_pdb[key]).translate_crystal(pdb=mix_pdb[key],translate=[0,0,4000])
            chimera_=chimera.Chimera(path_wd+'/'+mix_pdb[key])

            print('-- Generate '+str(key)+' system from '+str(mix_pdb[key])+' --')
            print('-- Please wait, this may take some time ... --')
            chimera_.matrixset(pdb=mix_pdb[key],crystalcontacts=system.crystalcontacts.crystalcontacts_file,
                    system_size=system_size,fibril_length=fibril_length)
        
            print('-- Cut '+str(key)+' system to '+str(fibril_length)+' nm --')
            system_=matrixset_system(system=system,crystalcontacts_file=system.crystalcontacts.crystalcontacts_file)
            
            print('-- Add caps --')
            subprocess.run('rm -r '+path_wd+'/'+key,shell=True)
            subprocess.run('mkdir '+path_wd+'/'+key,shell=True)
            cap_system(system=system,crosslink_type=key)

        print('-- Mix system --')
        mix_=mix.Mix(ratio_mix=mix_setup,system=system)
        system_=mix_.add_mix(system=system)

    elif ratio_mix==None and connect_file!=None:

        for pdb in pdb_files[1:]:
            crystal.Crystal(pdb=pdb).translate_crystal(pdb=pdb,translate=[0,0,4000])
        connect_file=connect_file.replace('.txt','')

        print('-- Mix system --')
        print('-- NOTE: Make sure each crosslink-type in '+str(connect_file)+' (D,T,DT,TD) is provided --')
        mix_=mix.Mix(system=system,connect_mix=connect_file)
        system_=mix_.get_mix_from_connect_file(connect_file=connect_file)
    else:
        print('Error: Either provide a connect-file with the crosslink-type (e.g. 1.caps.pdb 2.caps.pdb ; D) OR \n'+
              'provide the -ratio_mix and -files_mix to generate a differently crosslinked microfibril \n'+
              'in both cases -files_mix flag has to contain the pdb-file with the respective crosslink types '+
              'NOTE: Make sure that each crosslink type (D,T,DT,TD) from connect-file is provided!')
        
    connect.Connect(system=system_).write_connect(system=system_,connect_file=connect_file)

    system_.write_pdb(pdb_out=pdb_out,fibril_length=fibril_length)

    return system_

def build_system(crystal: crystal.Crystal,
                 crystalcontacts: crystalcontacts.CrystalContacts) -> system.System:
    """
    
    build a system of models, i.e., a microfibril of triple helices
    
    """
    system_=system.System(crystal=crystal,crystalcontacts=crystalcontacts)

    transformation=system_.crystalcontacts.read_t_matrix()
    unit_cell={ k:system_.crystal.get_s_matrix(t_matrix=transformation[k]) for k in transformation }

    for key_m in transformation:
        model_=model.Model(id=key_m,transformation=transformation[key_m],unit_cell=unit_cell[key_m],pdb_file=crystal.pdb_file)
        system_.add_model(model=model_)
    return system_

def connect_system(system: system.System,connect_file=str) -> tuple[system.System, connect.Connect]:
    """
    
    identify crosslink connections between models within the system 
    
    """
    connect_=connect.Connect(system=system,connect_file=connect_file)
    system_connect=connect_.run_connect(system=system)
    for key_m in system_connect:
        system.get_model(model_id=key_m).add_connect(connect_id=key_m,connect=system_connect[key_m])
    return system,connect_

def cap_system(system: system.System, crosslink_type: str) -> caps.Caps:
    """
    
    cap each model of system
    
    """
    caps_=caps.Caps(system=system)
    for idx in system.get_models():
        caps_.read_residues(pdb_id=int(idx))
        caps_.add_caps(pdb_id=int(idx),crosslink_type=crosslink_type)
    return caps_

def matrixset_system(system: system.System,crystalcontacts_file=str) -> system.System:
    """
    
    set system after cutting the fibril to its desired length.
    NOTE: Some models are deleted, since they are outside of the fibril length
    
    """
    contacts=[float(i.split(' ')[1]) for i in open(crystalcontacts_file.replace('_opt','')+'_id.txt','r').readlines()]
    for model in system.get_models():
        if model not in contacts: 
            system.delete_model(model_id=model)
        elif system.get_model(model_id=model).connect!=None:
            for connect in system.get_model(model_id=model).connect:
                if connect not in contacts: system.get_model(model_id=model).delete_connect(connect_id=connect)
    return system

def build_from_contactdistance(path_wd=str,pdb_file=None,contact_distance=float,crystalcontacts_file=str,solution_space=[],connect_file=str,
                               chimera=chimera.Chimera,crystal=crystal.Crystal) -> tuple[system.System, crystalcontacts.CrystalContacts, connect.Connect]:
    """

    generate system of models based on contact distance and PDB-file
    
    """
    path_pdb_file=path_wd+'/'+pdb_file

    crystalcontacts_file='crystalcontacts_from_colbuilder' # default name for crystalcontacts file
    connect_file='connect_from_colbuilder' # default name for connect file

    print('-- Get CrystalContacts for contact distance '+str(contact_distance)+' A --')
    print('-- Please wait, this may take some time ... --')
    chimera.matrixget(pdb=path_pdb_file,contact_distance=contact_distance,
                       crystalcontacts=crystalcontacts_file)

    print('-- Write '+str(crystalcontacts_file)+' --')
    crystalcontacts_=crystalcontacts.CrystalContacts(crystalcontacts_file)

    print('-- Build system --')
    system_=build_system(crystal=crystal,crystalcontacts=crystalcontacts_)
        
    print('-- Connect system --')
    system_,connect_=connect_system(system=system_,connect_file=connect_file)

    print('-- Optimize system --')
    system_=optimize.Optimizer(system=system_,solution_space=solution_space).run_optimize(system=system_,connect=connect_)
    system_,connect_=connect_system(system=system_,connect_file=connect_file)
    
    crystalcontacts_.crystalcontacts_file=crystalcontacts_file+'_opt'

    return system_,crystalcontacts_,connect_

def build_from_crystalcontacts(crystalcontacts_file=str,crystal=crystal.Crystal,solution_space=[],connect_file=str,
                               crystalcontacts_optimize=bool) -> tuple[system.System, crystalcontacts.CrystalContacts, connect.Connect]:
    """
    
    generate system of models based on CrystalContacts and PDB-file 

    """
    crystalcontacts_file=crystalcontacts_file.replace('.txt','')
    if connect_file==None: 
        bool_external_connect=False 
        connect_file='connect_from_colbuilder' # default name for connect file
    else:
        bool_external_connect=True
        connect_file=connect_file.replace('.txt','')

    crystalcontacts_=crystalcontacts.CrystalContacts(crystalcontacts_file)
    
    print('-- Build system --')
    system_=build_system(crystal=crystal,crystalcontacts=crystalcontacts_)

    print('-- Connect system --')
    system_,connect_=connect_system(system=system_,connect_file=connect_file)
   
    if crystalcontacts_optimize:
        print('-- Optimize system --')
        system_=optimize.Optimizer(system=system_,solution_space=solution_space).run_optimize(system=system_,connect=connect_)
        system_,connect_=connect_system(system=system_)

        crystalcontacts_.crystalcontacts_file=crystalcontacts_.crystalcontacts_file+'_opt'

    if bool_external_connect:
        system_=connect_.get_external_connect_file(system=system_,connect_file=connect_file)

    return system_,crystalcontacts_,connect_

def build_fibril(path_wd=str,pdb_file=str) -> system.System:
    """"
    
    build a system from colbuilder 1.0 fibril
    
    """
    system_=system.System(pdb_fibril=pdb_file)

    print('-- Read fibril '+str(pdb_file)+' from colbuilder 1.0 --')
    fibril_=fibril.Fibril(system=system_,pdb_file=pdb_file)

    print('-- Separate system --')
    fibril_.seperate_system(pdb_file=pdb_file)

    print('-- Build system --')
    system_=fibril_.build_system(system=system_)

    print('-- Connect system --')
    system_,_=connect_system(system=system_)
    
    print('-- Write connect_from_colbuilder.txt --')
    fibril_.write_connect(system=system_,connect_file='connect_from_colbuilder')

    return system_

def build_from_pdb(path_wd=str,pdb_file=None,contact_distance=float,crystalcontacts_file=str,
                    chimera=chimera.Chimera,crystal=crystal.Crystal) -> tuple[system.System, crystalcontacts.CrystalContacts, connect.Connect]:
    """

    generate system of models based on contact distance and PDB-file
    
    """
    path_pdb_file=path_wd+'/'+pdb_file

    crystalcontacts_file='crystalcontacts_from_colbuilder' # default name for crystalcontacts file

    print('-- Prepare topology for single PDB-file --')
    chimera.matrixget(pdb=path_pdb_file,contact_distance=contact_distance,
                       crystalcontacts=crystalcontacts_file)
    
    crystalcontacts_=crystalcontacts.CrystalContacts(crystalcontacts_file)

    print('-- Build system --')
    system_=build_system(crystal=crystal,crystalcontacts=crystalcontacts_)
    system_,connect_=connect_system(system=system_)
        
    return system_,crystalcontacts_,connect_
