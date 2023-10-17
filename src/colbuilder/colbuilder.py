
import argparse
from pathlib import Path

from colbuilder.geometry.main_geometry import build_geometry, mix_geometry, mutate_geometry, build_fibril
from colbuilder.topology.main_topology import build_topology
    
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
                        help='lengh of microfibril (default: 334 nm)',default=334)
    parser.add_argument('-contacts','--crystalcontacts_file', required=False, 
                        help='read crystalcontacts from file (default: crystalcontacts)',default='crystalcontacts')
    parser.add_argument('-connect','--connect_file', required=False, 
                        help='read external crystalcontacts-connect file  (default: crystalcontacts-file_connect)',default='')
    parser.add_argument('-optimize','--crystalcontacts_optimize', action='store_true', 
                        help='optimize crystalcontacts (default: False)',default=False)
    parser.add_argument('-geometry','--geometry_generator', action='store_true', 
                        help='generate geometry files (default: False)',default=False)
    parser.add_argument('-space','--solution_space', nargs='+', required=False,
                        help='solution space of optimisation problem [ d_x d_y d_z ] (default: [1 1 1] )',default=[1,1,1])
    
    parser.add_argument('-fibril', '--fibril', required=False, 
                        help='PDB-file of colbuilder 1 fibril',default=None)
    
    parser.add_argument('-mix','--setup_mix', required=False,nargs='+',
                        help=("ratio for mix-crosslink setup, e.g. 0.7 T; 0.3 D -> -mix t:70 d:30. Please use -f_mix flag to input pdb-files in the exact same order"),default=None)
    parser.add_argument('-fmix','--files_mix', required=False,nargs='+',
                        help=("PDB-files with different crosslink-types, e.g. 0.7 T; 0.3 D -> fmix Rat-T.pdb Rat-D.pdb. Please make sure that -fmix has the exact same order as mix-setup -mix."),default=[])
    parser.add_argument('-mutate','--setup_mutate', required=False,
                        help=("ratio of mutated crosslinks, e.g. -mutate 0.25 -> 0.25 mutated, values between 0 to 0.5"),default=None)
    
    parser.add_argument('-topology','--topology_generator', action='store_true', 
                        help='generate topology files (default: False)',default=False)
    parser.add_argument('-go','--go_eps', required=False,
                        help=("specifiy potential well of go-like potential (default: 9.414)"),default='9.414')
    parser.add_argument('-p','--topology_file', required=False,
                        help=("specifiy name of topology file (default: system.top)"),default='system.top')
    parser.add_argument('-ff','--force_field', required=False,
                        help=("specifiy force field to be used, e.g. -ff amber99 OR -ff martini3"),default=None)
        
    args=parser.parse_args()
    if args.connect_file=='': args.connect_file=str(args.crystalcontacts_file).replace('.txt','')+'_connect.txt'

    print('-- Colbuilder 2.0 --')

    if args.file==None and args.files_mix!=[]: args.file=args.files_mix[0]

    # Build Geometry of Microfibril
    if args.fibril==None:
        system_=build_geometry(path_wd=str(args.working_directory),
                        pdb_file=str(args.file).replace('.pdb',''),
                        contact_distance=args.contact_distance,
                        crystalcontacts_file=str(args.crystalcontacts_file).replace('.txt',''),
                        crystalcontacts_optimize=args.crystalcontacts_optimize,
                        connect_file=str(args.connect_file).replace('.txt',''),
                        solution_space=args.solution_space,
                        fibril_length=float(args.fibril_length),
                        geometry=args.geometry_generator,
                        pdb_out=args.output)
    
    if args.file==None and args.fibril!=None:
        system_=build_fibril(path_wd=str(args.working_directory),
                            pdb_file=args.fibril)

    # Mix-System
    if args.setup_mix!=None: 
        system_=mix_geometry(path_wd=str(args.working_directory),
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

    # Build Topology for System
    if args.topology_generator==True:
        system_=build_topology(system=system_,
                           force_field=args.force_field,
                           top_file=args.output+'.top',
                           gro_file=args.output+'.gro',
                           go_epsilon=args.go_eps)

if __name__ == '__main__':
    colbuilder()