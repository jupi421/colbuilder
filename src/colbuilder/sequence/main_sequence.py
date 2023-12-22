from colbuilder.sequence import alignment,model


def build_sequence(path_wd=str,pdb_file=None,collagen_type=int,
                      dict_chain={},register=[],crosslink={}):
    """
    
    build sequence from an uncrossed collagen triple helix
    
    """
    print('-- Building the sequence of the Collagen Triple Helix --')

    print('-- Perform the sequence alignment step --')
    model_=model.Model(pdb_filename=pdb_file)
    model_.read_pdb(pdb_filename=pdb_file)
    print(model_.pdb['atom_id'])

