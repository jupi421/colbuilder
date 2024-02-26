import argparse
from pathlib import Path

from colbuilder.geometry.main_geometry import build_geometry, mix_geometry, mutate_geometry, build_fibril
from colbuilder.topology.main_topology import build_topology
from colbuilder.sequence.main_sequence import build_sequence
    
def colbuilder():
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-f', '--file', required=False, 
                        help='PDB-input file for single triple helix or colbuilder 1.0 fibril',default=None)
    parser.add_argument('-o', '--output', required=False, 
                        help='Name for PDB-file of microfibril ',default='collagen_fibril')
    parser.add_argument('-wd','--working_directory', required=False, 
                        help='set working directory ',default=Path.cwd())
    parser.add_argument('-dc','--contact_distance', required=False, 
                        help='contact distance as input for radial size of microfibril, e.g. 10 to 60 ',default=None)
    parser.add_argument('-length','--fibril_length', required=False, 
                        help='lengh of microfibril ',default=334)
    parser.add_argument('-contacts','--crystalcontacts_file', required=False, 
                        help='read crystalcontacts from file ',default='crystalcontacts_from_colbuilder')
    parser.add_argument('-connect','--connect_file', required=False, 
                        help='read connect between contacts from file',default='connect_from_colbuilder')
    parser.add_argument('-optimize','--crystalcontacts_optimize', action='store_true', 
                        help='optimize crystalcontacts ',default=False)
    parser.add_argument('-geometry','--geometry_generator', action='store_true', 
                        help='generate geometry files ',default=False)
    parser.add_argument('-space','--solution_space', nargs='+', required=False,
                        help='solution space of optimisation problem [ d_x d_y d_z ] ',default=[1,1,1])
    parser.add_argument('-fibril', '--fibril', required=False, action='store_true', 
                        help='Bool argument to generate topology for colbuilder 1.0 67nm-long fibril ',default=False)
    
    parser.add_argument('-ratio_mix','--ratio_mix', required=False,nargs='+',
                        help=("ratio for mix-crosslink setup, e.g. 0.7 T; 0.3 D -> -mix T:70 D:30. Please use -files_mix flag to input pdb-files in the exact same order"),default=None)
    parser.add_argument('-files_mix','--files_mix', required=False,nargs='+',
                        help=("PDB-files with different crosslink-types, e.g. 0.7 T; 0.3 D -> fmix Rat-T.pdb Rat-D.pdb."+ 
                              "If the ratio_mix is provided, please make sure that files_mix has the exact same order as -ratio_mix flag OR"+
                               "If connect_mix information is provided make sure to provide each triple helix crosslink type as input for -files_mix."),default=None)
    parser.add_argument('-connect_mix','--connect_mix', required=False,
                         help=("Provide connect file with mixture of triple helices within microfibril"),default=None)
    
    parser.add_argument('-mutate','--setup_mutate', required=False,
                        help=("ratio of mutated crosslinks, e.g. -mutate 0.25 -> 0.25 mutated, values between 0 to 0.5"),default=None)
    
    parser.add_argument('-topology','--topology_generator', action='store_true', 
                        help='generate topology files ',default=False)
    parser.add_argument('-go','--go_eps', required=False,
                        help=("specifiy potential well of go-like potential "),default='9.414')
    parser.add_argument('-p','--topology_file', required=False,
                        help=("specifiy name of topology file "),default='system.top')
    parser.add_argument('-ff','--force_field', required=False,
                        help=("specifiy force field to be used, e.g. -ff amber99 OR -ff martini3"),default=None)
    
    parser.add_argument('-sequence','--sequence_generator', action='store_true', 
                        help='generate triple helix from sequence ',default=False)
    parser.add_argument('-type','--collagen_type', required=False,
                        help=("specifiy type of collagen molecule "),default=1)
    parser.add_argument('-crosslink','--crosslink_topology', required=False,
                        help=("specifiy crosslink types of the triple helix "),default=['no','no'])
    parser.add_argument('-register','--register_topology', required=False,
                        help=("specifiy register of chains, e.g., A,A,C "),default='A,B,C')
    parser.add_argument('-chain','--chain_id', required=False,
                        help=("specifiy the chain id: chain 1=A, chain 2=C, chain 3=B "),default='A,B,C')
    parser.add_argument('-ensemble','--ensemble', required=False,
                        help=("Ensemble of modeller runs "),default=3)
        
    args=parser.parse_args()
    if args.connect_file=='': args.connect_file=str(args.crystalcontacts_file).replace('.txt','_connect.txt')

    print('-- Colbuilder 2.0 --')

    if args.file==None and args.files_mix!=[]: args.file=args.files_mix[0]

    # Build Triple Helix from sequence
    if args.sequence_generator==True:
        system_=build_sequence(path_wd=str(args.working_directory),
                        pdb_file=str(args.file).replace('.pdb',''),
                        collagen_type=args.collagen_type,
                        dict_chain={ i: args.chain_id[i] for i in range(len(args.chain_id.split(',')))},
                        register=[str(i) for i in args.register_topology.split(',')],
                        crosslink=args.crosslink_topology,
                        ensemble=args.ensemble)

    # Build Geometry of Microfibril
    if args.fibril==False:
        system_=build_geometry(path_wd=str(args.working_directory),
                        pdb_file=str(args.file).replace('.pdb',''),
                        contact_distance=args.contact_distance,
                        crystalcontacts_file=str(args.crystalcontacts_file).replace('.txt',''),
                        crystalcontacts_optimize=args.crystalcontacts_optimize,
                        connect_file=str(args.connect_file).replace('.txt',''),
                        solution_space=args.solution_space,
                        fibril_length=float(args.fibril_length),
                        geometry=args.geometry_generator,
                        pdb_out=str(args.output).replace('.pdb',''))
    
    if args.fibril==True:
        system_=build_fibril(path_wd=str(args.working_directory),
                            pdb_file=args.file,
                            connect_file=str(args.connect_file).replace('.txt',''))

    # Mix-System
    if args.files_mix!=None: 
        system_=mix_geometry(path_wd=str(args.working_directory),
                            fibril_length=float(args.fibril_length),
                            pdb_files=[str(file).replace('.pdb','') for file in args.files_mix],
                            ratio_mix=args.ratio_mix,
                            connect_file_mix=args.connect_mix.replace('.txt',''),
                            system=system_,
                            pdb_out=str(args.output).replace('.pdb','_mix'))
        
    # Mutate-System
    if args.setup_mutate!=None:
        system_=mutate_geometry(path_wd=str(args.working_directory),
                                setup_mutate=args.setup_mutate,
                                system=system_,
                                fibril_length=float(args.fibril_length),
                                pdb_out=str(args.output).replace('.pdb','_mut'))

    # Build Topology for System
    if args.topology_generator==True:
        system_=build_topology(system=system_,
                           force_field=args.force_field,
                           top_file=args.output+'.top',
                           gro_file=str(args.output).replace('.pdb','.gro'),
                           go_epsilon=args.go_eps)

if __name__ == '__main__':
    colbuilder()