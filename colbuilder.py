#!/hits/fast/mbm/broszms/software/conda/envs/colbuilder/bin/python3.9

import os
from setuptools import setup
#
import argparse
from colbuilder.geometry.build_geometry import build_geometry

def colbuilder():
    
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

    # Build Geometry of Microfibril
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
    colbuilder()