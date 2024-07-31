import pytest
import asyncio
from pathlib import Path
import pandas as pd
import tempfile
import shutil
import os
from time import sleep

from colbuilder.sequence.main_sequence import build_sequence
from colbuilder.sequence.alignment import Alignment, align_sequences
from colbuilder.sequence.modeller import ModellerWrapper, run_modeller
from colbuilder.sequence.mutate_crosslinks import apply_crosslinks, parse_crosslink_info, rename_residue_in_pdb
from colbuilder.core.utils.config import ColbuilderConfig

# Mock data and fixtures
@pytest.fixture
def mock_config():
    return ColbuilderConfig(
        fasta_file="test.fasta",
        TEMPLATE_FASTA_PATH="template.fasta",
        TEMPLATE_PDB_PATH="template.pdb",
        RESTYP_LIB_PATH="restyp.lib",
        TOP_HEAV_LIB_PATH="top_heav.lib",
        PAR_MOD_LIB_PATH="par_mod.lib",
        species="test_species",
        n_term_type="test_n",
        c_term_type="test_c",
        n_term_combination="test_n_combo",
        c_term_combination="test_c_combo",
        CROSSLINKS_FILE="crosslinks.csv",
        crosslink=True,
        debug=True,
        pdb_first_line="TEST PDB FIRST LINE"
    )

@pytest.fixture
def temp_dir():
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)

# Test main_sequence.py
@pytest.mark.asyncio
async def test_build_sequence(mock_config, temp_dir, mocker):
    mock_config.fasta_file = Path(temp_dir) / "test.fasta"
    mock_config.CROSSLINKS_FILE = Path(temp_dir) / "crosslinks.csv"
    mock_config.species = "test_species"
    mock_config.n_term_type = "test_n"
    mock_config.c_term_type = "test_c"
    mock_config.n_term_combination = "test_n_combo"
    mock_config.c_term_combination = "test_c_combo"
    mock_config.crosslink = True
    mock_config.debug = True
    mock_config.pdb_first_line = "TEST PDB FIRST LINE"

    # Create mock input files
    with open(mock_config.fasta_file, "w") as f:
        f.write(">test\nABCDEF\n")

    with open(mock_config.CROSSLINKS_FILE, "w") as f:
        f.write("species,terminal,type,combination,R1,P1,R2,P2,R3,P3\n")
        f.write("test_species,N,test_n,test_n_combo,DHLNL,87.A,NONE,,NONE,\n")

    # Mock the debug_output directory
    debug_output = Path.cwd() / 'debug_output'
    debug_output.mkdir(exist_ok=True)

    # Create mock PDB and alignment files
    mock_output_pdb = debug_output / "test_final_model.pdb"
    mock_crosslinked_pdb = debug_output / "test_N_test_n_C_NONE_temp.pdb"
    mock_formatted_pdb = debug_output / "test_N_test_n_C_NONE.pdb"
    mock_alignment_fasta = debug_output / "test_alignment.fasta"
    mock_modeller_ali = debug_output / "test_modeller.ali"

    for pdb_file in [mock_output_pdb, mock_crosslinked_pdb, mock_formatted_pdb]:
        with open(pdb_file, "w") as f:
            f.write("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n")

    with open(mock_alignment_fasta, "w") as f:
        f.write(">test_aligned\nABCDEF\n")

    with open(mock_modeller_ali, "w") as f:
        f.write("Mock MODELLER alignment file content")

    # Mock functions
    mocker.patch("colbuilder.sequence.main_sequence.align_sequences", return_value=(str(mock_alignment_fasta), str(mock_modeller_ali)))
    mocker.patch("colbuilder.sequence.main_sequence.run_modeller", return_value=str(mock_output_pdb))
    mocker.patch("colbuilder.sequence.main_sequence.apply_crosslinks", return_value=str(mock_crosslinked_pdb))
    mocker.patch("colbuilder.sequence.main_sequence.format_pdb", return_value=None)

    # Run the function
    msa_output, final_output = await build_sequence(mock_config)

    # Assertions
    assert isinstance(msa_output, Path)
    assert isinstance(final_output, Path)
    assert final_output.exists()
    assert final_output.name == "test_N_test_n_C_NONE.pdb"

    # Check if the formatted PDB file was created
    assert mock_formatted_pdb.exists()

    # Verify the content of the formatted PDB file
    with open(mock_formatted_pdb, "r") as f:
        content = f.read()
        assert "ATOM" in content

    # Check if the MSA file was created and copied
    assert mock_alignment_fasta.exists()
    assert msa_output.exists()
    assert msa_output.name == "test_alignment.fasta"

    # Clean up
    shutil.rmtree(debug_output)

# Test alignment.py
def test_alignment_class(temp_dir, mock_config):
    input_fasta = Path(temp_dir) / "input.fasta"
    template_fasta = Path(temp_dir) / "template.fasta"
    
    with open(input_fasta, "w") as f:
        f.write(">input\nABCDEF\n")
    with open(template_fasta, "w") as f:
        f.write(">template\nXYZABC\n")
    
    alignment = Alignment(input_fasta, template_fasta, "test_prefix", Path(temp_dir) / "template.pdb")
    
    assert len(alignment.input_sequences) == 1
    assert len(alignment.template_sequences) == 1

def test_align_sequences(temp_dir, mock_config, mocker):
    input_fasta = Path(temp_dir) / "input.fasta"
    template_fasta = Path(temp_dir) / "template.fasta"
    template_pdb = Path(temp_dir) / "template.pdb"
    
    with open(input_fasta, "w") as f:
        f.write(">input\nABCDEF\n")
    with open(template_fasta, "w") as f:
        f.write(">template\nXYZABC\n")
    with open(template_pdb, "w") as f:
        f.write("ATOM      1  N   ALA A   1      -0.525   1.362   0.000  1.00  0.00           N  \n")

    mocker.patch("colbuilder.sequence.alignment.Alignment.align_sequences", return_value=("msa.fasta", "modeller.ali"))
    
    result = align_sequences(input_fasta, template_fasta, "test_prefix", template_pdb)
    
    assert len(result) == 2
    assert isinstance(result[0], str)
    assert isinstance(result[1], str)

# Test modeller.py
def test_modeller_wrapper(temp_dir, mock_config, mocker):
    mocker.patch("colbuilder.sequence.modeller.Environ")
    mocker.patch("colbuilder.sequence.modeller.AutoModel")
    mocker.patch("os.rename")
    
    wrapper = ModellerWrapper(
        str(Path(temp_dir) / "aligned.ali"),
        str(Path(temp_dir) / "template.pdb"),
        "test_prefix",
        str(Path(temp_dir) / "restyp.lib"),
        str(Path(temp_dir) / "top_heav.lib"),
        str(Path(temp_dir) / "par_mod.lib")
    )
    
    wrapper.run_modeller()
    
    assert wrapper.output_pdb is not None

def test_run_modeller(temp_dir, mock_config, mocker):
    mocker.patch("colbuilder.sequence.modeller.ModellerWrapper.execute_modeller")
    mocker.patch("colbuilder.sequence.modeller.ModellerWrapper.output_pdb", new_callable=mocker.PropertyMock, return_value="output.pdb")
    
    result = run_modeller(
        str(Path(temp_dir) / "aligned.ali"),
        str(Path(temp_dir) / "template.pdb"),
        "test_prefix",
        str(Path(temp_dir) / "restyp.lib"),
        str(Path(temp_dir) / "top_heav.lib"),
        str(Path(temp_dir) / "par_mod.lib")
    )
    
    assert result == "output.pdb"

# Test mutate_crosslinks.py
def test_rename_residue_in_pdb(temp_dir):
    pdb_content = "ATOM      1  N   ALA A   1      -0.525   1.362   0.000  1.00  0.00           N  \n"
    pdb_file = Path(temp_dir) / "test.pdb"
    with open(pdb_file, "w") as f:
        f.write(pdb_content)
    
    rename_residue_in_pdb(str(pdb_file), "A", 1, "GLY")
    
    with open(pdb_file, "r") as f:
        new_content = f.read()
    
    assert "GLY A   1" in new_content

def test_parse_crosslink_info():
    crosslink_row = pd.Series({
        "R1": "DHLNL", "P1": "87.A",
        "R2": "NONE", "P2": None,
        "R3": "NONE", "P3": None
    })
    
    result = parse_crosslink_info(crosslink_row)
    
    assert len(result) == 1
    assert result[0] == ("DHLNL", "A", 87)

def test_apply_crosslinks(temp_dir, mock_config, mocker):
    mocker.patch("colbuilder.sequence.mutate_crosslinks.Environ")
    mocker.patch("colbuilder.sequence.mutate_crosslinks.complete_pdb")
    mocker.patch("colbuilder.sequence.mutate_crosslinks.rename_residue_in_pdb")
    
    input_pdb = Path(temp_dir) / "input.pdb"
    output_pdb = Path(temp_dir) / "output.pdb"
    with open(input_pdb, "w") as f:
        f.write("ATOM      1  N   ALA A   1      -0.525   1.362   0.000  1.00  0.00           N  \n")
    
    n_crosslink = pd.Series({"R1": "DHLNL", "P1": "87.A"})
    c_crosslink = None
    
    result = apply_crosslinks(str(input_pdb), str(output_pdb), n_crosslink, c_crosslink, mock_config)
    
    assert result == str(output_pdb)

if __name__ == "__main__":
    pytest.main([__file__])