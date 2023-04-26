#!/hits/fast/mbm/broszms/software/conda/envs/colbuilder/bin/python3.9

import os
import sys
from setuptools import setup
import typing
#
import argparse
import logging
from geometry import (
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

        print('-- Save CrystalContacts of Systems in '+str(crystalcontacts_file)+'.txt --')
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

    elif pdb_file!=None and contact_distance==None and crystalcontacts_file!='chimera_crystalcontacts':
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

    print('-- Microfibril was generated, now N-, and C- termini are capped ... --')
    caps_=cap_system(system=system_)


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
    caps_=caps.Caps(system_size=system.size_models)
    for model_id in range(system.size_models):
        caps_.add_caps(model_id=model_id)

    return caps_
    print('hhelo')
   # for model in system.size_system


def generate_atomistic_topology(system : system.System,crystalcontacts_file=str):
    """
    
    Generates atomistis topology based on connect information of system
    
    """


    print('hello world')




    return print('Atomistic topology written')



def main():
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    print('-- Colbuilder --')
    print('-- Read input parameters --')
    
    # Input arguments 
    parser.add_argument('-f', '--file', required=True, 
                        help='PDB-file of single triple helix (PDB)')
    parser.add_argument('-wd','--working_directory', required=False, 
                        help='set working directory',default=os.getcwd())
    parser.add_argument('-dc','--contact_distance', required=False, 
                        help='contact distance as input for crystalcontacts command',default=None)
    parser.add_argument('-cut-off','--cut_off', required=False, 
                        help='cut-off Microfibril 300 +/- 15 [nm] ',default=315)
    parser.add_argument('-contacts','--crystalcontacts_file', required=False, 
                        help='read user-specific crystalcontacts information from *.txt file',default='chimera_crystalcontacts')
    parser.add_argument('-optimize','--crystalcontacts_optimize', required=False, 
                        help='optimize user-specified crystalcontacts information',default=False)
    args=parser.parse_args()

    # TODO: PATH LOGIC    
    system_=build_geometry(path_wd=str(args.working_directory),
                            pdb_file=str(args.file).replace('.pdb',''),
                            contact_distance=args.contact_distance,
                            crystalcontacts_file=str(args.crystalcontacts_file).replace('.txt',''),
                            crystalcontacts_optimize=args.crystalcontacts_optimize,
                            cut_off=float(args.cut_off))

    # TODO: AA FF
   # system_aa_=generate_atomistic_topology(system=system_,crystalcontacts_file=str(args.crystalcontacts_file).replace('.txt',''))
    
    # TODO CG FF
    
    # MAP GO-Model
#    map_go_itp.run_map_go_itp(path,file_name)


if __name__ == '__main__':
    main()
