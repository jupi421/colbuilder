import argparse
from pathlib import Path

from colbuilder.geometry.main_geometry import build_geometry, mix_geometry, replace_geometry, build_fibril
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
                        help='read crystalcontacts from file ',default=None)
    parser.add_argument('-connect','--connect_file', required=False, 
                        help='read connect between contacts from file',default=None)
    parser.add_argument('-optimize','--crystalcontacts_optimize', action='store_true', 
                        help='optimize crystalcontacts ',default=False)
    parser.add_argument('-geometry','--geometry_generator', action='store_true', 
                        help='generate geometry files ',default=False)
    parser.add_argument('-space','--solution_space', nargs='+', required=False,
                        help='solution space of optimisation problem [ d_x d_y d_z ] ',default=[1,1,1])
    parser.add_argument('-fibril', '--fibril', required=False, action='store_true', 
                        help='Bool argument to generate topology for colbuilder 1.0 67nm-long fibril ',default=False)
    
    parser.add_argument('-mix','--mix_bool', required=False,action='store_true',
                         help=("Set -mix flag to generate a mixed crosslinked microfibril"),default=False)
    parser.add_argument('-ratio_mix','--ratio_mix', required=False,nargs='+',
                        help=("Ratio for mix-crosslink setup: -ratio_mix T:70 D:30\n"+
                               "Provide files at -files_mix flag in same order as for ratio_mix"),default=None)
    parser.add_argument('-files_mix','--files_mix', required=False,nargs='+',
                        help=("PDB-files with different crosslink-types: -files_mix Rat-T.pdb Rat-D.pdb\n"+ 
                              "If the ratio_mix is provided, make sure that -files_mix has the same order as -ratio_mix OR\n"+
                              "If connect information is provided, provide each type of crosslinked triple helix as input for -files_mix."),default=[])
    
    parser.add_argument('-replace','--replace_bool', required=False,action='store_true',
                         help=("Set -replace flag to generate a microfibril with less crosslinks"),default=False)
    parser.add_argument('-ratio_replace','--ratio_replace', required=False,
                        help=("Ratio of crosslinks to be replaced with Lysines: -ratio_replace 25 means that 25"+
                              " crosslinks are replaced with Lysines (range: 0 to 50"+")"),default=None)
    parser.add_argument('-replace_file','--replace_file', required=False,
                        help=("File with information about crosslinks to be replaced with Lysine. Each crosslink to be replaces should be noted according to this syntax: "+
                              " molecule_id, residue_name, residue_id, chain_id (e.g. 1.caps.pdb LY2 1046 A)"),default=None)
    
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

    if args.mix_bool==True and args.files_mix==[]: args.files_mix=[args.file]
    if args.mix_bool==False and args.files_mix!=[]: print("Error: Please set -mix flag to obtain a mixed structure or topology."); exit()

    # Build a triple helix from amino acid sequence
    if args.sequence_generator:
    final_pdb = build_sequence(path_wd=args.working_directory,
                               fasta_file=args.fasta_file,
                               dict_chain={ i: args.chain_id[i] for i in range(len(args.chain_id.split(',')))},
                               crosslink=args.crosslink_topology)
    
    print(f"Final PDB file created: {final_pdb}")

    # Build a system of models, i.e., the geometry, for a long microfibril
    if args.fibril==False and args.sequence_generator==False:
        system_=build_geometry(path_wd=str(args.working_directory),
                        pdb_file=str(args.file).replace('.pdb',''),
                        contact_distance=args.contact_distance,
                        crystalcontacts_file=args.crystalcontacts_file,
                        crystalcontacts_optimize=args.crystalcontacts_optimize,
                        connect_file=args.connect_file,
                        solution_space=args.solution_space,
                        fibril_length=float(args.fibril_length),
                        geometry=args.geometry_generator,
                        pdb_out=str(args.output).replace('.pdb',''))
    
    # Build a system of models for the 67-nm long collagen D-Band from colbuilder1
    if args.fibril==True and args.sequence_generator==False:
        system_=build_fibril(path_wd=str(args.working_directory),
                            pdb_file=args.file,
                            connect_file=args.connect_file)

    # Mix divalent and trivalent crosslinks within system to alter crosslink density
    if args.mix_bool==True and args.replace_bool==False: 
        system_=mix_geometry(path_wd=str(args.working_directory),
                            fibril_length=float(args.fibril_length),
                            pdb_files=[str(file).replace('.pdb','') for file in args.files_mix],
                            ratio_mix=args.ratio_mix,
                            connect_file=args.connect_file,
                            system=system_,
                            pdb_out=str(args.output).replace('.pdb',''))

    # Replace crosslinks within system to reduce crosslink density
    elif args.replace_bool==True and args.mix_bool==False:
        system_=replace_geometry(path_wd=str(args.working_directory),
                                ratio_replace=args.ratio_replace,
                                system=system_,
                                fibril_length=float(args.fibril_length),
                                pdb_out=str(args.output).replace('.pdb',''),
                                replace_file=str(args.replace_file).replace('.txt',''))
    
    elif args.replace_bool==True and args.mix_bool==True:
        print('Error: -mix and -replace can not be used together, either mix crosslinks or replace them')

    # Build Topology for System
    if args.topology_generator==True and args.sequence_generator==False:
        system_=build_topology(system=system_,
                           force_field=args.force_field,
                           top_file=args.output+'.top',
                           gro_file=str(args.output).replace('.pdb','.gro'),
                           go_epsilon=args.go_eps)

if __name__ == '__main__':
    colbuilder()
