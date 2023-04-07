#!/hits/fast/mbm/broszms/software/conda/envs/colbuilder/bin/python3.10

# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 10:02:52 2023

@author: broszms
"""
import os
import sys
from setuptools import setup
#
# TODO: You can do better than this
path_colbuilder='/hits/fast/mbm/broszms/Collagen/colbuilder'
if path_colbuilder not in sys.path:
    sys.path.insert(0,'/hits/fast/mbm/broszms/Collagen/colbuilder')
#
import argparse
import logging
import time
from geometry import gen_coord
#from topology import map_go_itp
import analysis.triplehelix
#

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
    #
    # TODO: Better PATH logic -> path package
    path=str(args.working_directory)
    #
    # CrystalContacs are set -> Generate geometry
    if args.contact_distance!=None and args.file!=None:
        # Generate coordinate file for microfibril
        out=gen_coord.run_gen_coord(
            path,str(args.file).replace('.pdb',''),
            int(args.contact_distance),int(args.cut_off_fibril))
        #gen_coord.run_gen_coord(path,file_name,crystal_contacts,cut_off)
    #
    # TODO: AA FF
    #
    # TODO CG FF
    #
    # MAP GO-Model
#    map_go_itp.run_map_go_itp(path,file_name)

main()

if __name__ == '__main__':
    main()
    