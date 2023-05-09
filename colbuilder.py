
import os
from setuptools import setup
#
import argparse
from pathlib import Path
from colbuilder.geometry.main_geometry import build_geometry, mix_geometry, mutate_geometry

    
def colbuilder():
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-f', '--file', required=False, 
                        help='PDB-file of single triple helix',default=None)
    parser.add_argument('-o', '--output', required=False, 
                        help='Name for PDB-file of microfibril (default: collagen_fibril)',default='collagen_fibril')
    parser.add_argument('-wd','--working_directory', required=False, 
                        help='set working directory (default: cwd)',default=Path.cwd())
    parser.add_argument('-dc','--contact_distance', required=False, 
                        help='contact distance as input for radial size of microfibril, e.g. 10 to 60 (default: None)',default=None)
    parser.add_argument('-length','--fibril_length', required=False, 
                        help='lengh of microfibril (default: 315 nm)',default=315)
    parser.add_argument('-contacts','--crystalcontacts_file', required=False, 
                        help='read crystalcontacts from file (default: crystalcontacts)',default='crystalcontacts')
    parser.add_argument('-optimize','--crystalcontacts_optimize', action='store_true', 
                        help='optimize crystalcontacts (default: False)',default=False)
    parser.add_argument('-mix','--setup_mix', required=False,nargs='+',
                        help=("""ratio for mix-crosslink setup, e.g. 70% T; 30% D -> -mix T:70 D:30.
                              Please use -f_mix flag to input pdb-files in the exact same order"""),default=None)
    parser.add_argument('-fmix','--files_mix', required=False,nargs='+',
                        help=("""PDB-files with different crosslink-types, e.g. 70% T; 30 -> fmix Rat-T.pdb Rat-D.pdb
                        Please make sure that -fmix has the exact same order as mix-setup -mix."""),default=[])
    parser.add_argument('-mutate','--setup_mutate', required=False,
                        help=("ratio of mutated crosslinks, e.g. -mutate 25 means 25% mutated, values between 0 to 50%"),default=None)
    args=parser.parse_args()

    print('-- Colbuilder --')
    print('-- Read input parameters --')

    if args.file==None and args.files_mix!=[]: args.file=args.files_mix[0]

    # Build Geometry of Microfibril
    system_=build_geometry(path_wd=str(args.working_directory),
                        pdb_file=str(args.file).replace('.pdb',''),
                        contact_distance=args.contact_distance,
                        crystalcontacts_file=str(args.crystalcontacts_file).replace('.txt',''),
                        crystalcontacts_optimize=args.crystalcontacts_optimize,
                        fibril_length=float(args.fibril_length),
                        pdb_out=args.output)
    
    # Mix-System
    if args.setup_mix!=None: 
        system_=mix_geometry(path_wd=str(args.working_directory),
                            crystalcontacts_file=str(args.crystalcontacts_file).replace('.txt',''),
                            crystalcontacts_optimize=args.crystalcontacts_optimize,
                            fibril_length=float(args.fibril_length),
                            pdb_files=[str(file).replace('.pdb','') for file in args.files_mix],
                            setup_mix=args.setup_mix,
                            system=system_,
                            pdb_out=args.output+'_mix')
        
    # Mutate-System
    if args.setup_mutate!=None:
        system_=mutate_geometry(path_wd=str(args.working_directory),
                                setup_mutate=args.setup_mutate,
                                system=system_,
                                fibril_length=float(args.fibril_length),
                                pdb_out=args.output+'_mut')
    

    

    # TODO: AA FF
   # system_aa_=generate_atomistic_topology(system=system_,crystalcontacts_file=str(args.crystalcontacts_file).replace('.txt',''))
    
    # TODO CG FF
    
    # MAP GO-Model
#    map_go_itp.run_map_go_itp(path,file_name)

if __name__ == '__main__':
    colbuilder()
