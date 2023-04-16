#!/store/conda/envs/colbuilder/bin/python3.9

# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 10:02:52 2023

@author: broszms
"""
import os
import sys
from setuptools import setup
#
import argparse
import logging
from colbuilder.geometry.gen_coord import Fibril,Crystal
from colbuilder.geometry.get_crosslink import Crosslink
#
def main():
    #
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #
    # Input arguments 
    parser.add_argument('-f', '--file', required=False, 
                        help='PDB-file of single triple helix (PDB)')
    parser.add_argument('-wd','--working_directory', required=False, 
                        help='Set working directory',default=os.getcwd())
    parser.add_argument('-dc','--contact_distance', required=False, 
                        help='Number of Crystal Contacts (INT)',default=60)
    parser.add_argument('-cut-off','--cut_off_fibril', required=False, 
                        help='Cut-off Microfibril 300 +/- 15 [nm] ',default=315)
    args=parser.parse_args()
    
    # TODO: Better PATH logic -> path package
    path=str(args.working_directory)
    path_geo=__file__.replace('gen_coord.py','')
    pdb_file=str(args.file).replace('.pdb','')

    # CrystalContacs are set -> Generate geometry
    if args.contact_distance!=None and args.file!=None:
        # set fibril and get transformation & shift matrix from contacts
        fibril=Fibril(path_geo,path,pdb_file,int(args.contact_distance),int(args.cut_off_fibril))
        t_model,s_model=Crystal(fibril.contacts,fibril.pdb).run_system()
        # Translate crosslinks and search connection
        crosslink=Crosslink(fibril.pdb,t_model).run_system()
        

    # TODO: AA FF
    
    # TODO CG FF
    
    # MAP GO-Model
#    map_go_itp.run_map_go_itp(path,file_name)

main()

if __name__ == '__main__':
    main()
    
