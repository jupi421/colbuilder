import os
from typing import Any, Dict, Optional, Tuple, List
from pathlib import Path
from pydantic import BaseModel, Field
from enum import Enum
import yaml

class OperationMode(str, Enum):
    SEQUENCE = "sequence"
    GEOMETRY = "geometry"
    TOPOLOGY = "topology"
    FIBRIL = "fibril"
    MIX = "mix"
    REPLACE = "replace"
    BOTH = "both"

class ColbuilderConfig(BaseModel):
    mode: Optional[OperationMode] = Field(None, description="Operation mode")
    config_file: Optional[Path] = Field(None, description="YAML configuration file")
    sequence_generator: bool = Field(default=False, description="Run sequence generation")
    geometry_generator: bool = Field(default=False, description="Run geometry generation")
    fasta_file: Optional[str] = Field(None, description="FASTA file")
    output_prefix: Optional[str] = Field(None, description="Output prefix for sequence generation")
    species: Optional[str] = Field(None, description="Species name")
    crosslink: bool = Field(default=False, description="Apply crosslinks")
    n_term_type: Optional[str] = Field(None, description="N-terminal type")
    c_term_type: Optional[str] = Field(None, description="C-terminal type")
    n_term_combination: Optional[str] = Field(None, description="N-terminal combination")
    c_term_combination: Optional[str] = Field(None, description="C-terminal combination")
    file: Optional[Path] = Field(None, description="Input PDB file")
    output: Path = Field(default=Path("collagen_fibril"), description="Output file name")
    working_directory: Path = Field(default=Path.cwd(), description="Working directory")
    contact_distance: Optional[float] = Field(None, description="Contact distance for microfibril")
    fibril_length: float = Field(default=334, description="Length of microfibril")
    crystalcontacts_file: Optional[Path] = Field(None, description="Crystal contacts file")
    connect_file: Optional[Path] = Field(None, description="Connect file")
    crystalcontacts_optimize: bool = Field(default=False, description="Optimize crystal contacts")
    solution_space: Tuple[float, float, float] = Field(default=(1, 1, 1), description="Solution space")
    fibril: bool = Field(default=False, description="Generate topology for colbuilder 1.0 67nm-long fibril")
    mix_bool: bool = Field(default=False, description="Generate a mixed crosslinked microfibril")
    ratio_mix: Dict[str, int] = Field(default={}, description="Ratio for mix-crosslink setup")
    files_mix: Tuple[Path, ...] = Field(default=(), description="PDB files with different crosslink types")
    replace_bool: bool = Field(default=False, description="Generate a microfibril with less crosslinks")
    ratio_replace: Optional[float] = Field(None, description="Ratio of crosslinks to be replaced")
    replace_file: Optional[Path] = Field(None, description="File with crosslinks to be replaced")
    topology_generator: bool = Field(default=False, description="Generate topology files")
    go_eps: float = Field(default=9.414, description="Potential well of go-like potential")
    topology_file: Path = Field(default=Path("system.top"), description="Topology file name")
    force_field: Optional[str] = Field(None, description="Force field to be used")
    debug: bool = Field(default=False, description="Enable debug logging")
    # System settings ~constants
    pdb_first_line: Optional[str] = Field(default="CRYST1   39.970   26.950  677.900  89.24  94.59 105.58 P 1           2\n", description="Crystal contacts information")

    # Paths from your original config
    PROJECT_ROOT: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent)
    DATA_DIR: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data')
    HOMOLOGY_LIB_DIR: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'homology')
    TEMPLATE_PDB_PATH: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'homology' / "template.pdb")
    TEMPLATE_FASTA_PATH: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'homology' / "template.fasta")
    RESTYP_LIB_PATH: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'homology' / "modeller" / "restyp_mod.lib")
    TOP_HEAV_LIB_PATH: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'homology' / "modeller" / "top_heav_mod.lib")
    PAR_MOD_LIB_PATH: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'homology' / "modeller" / "par_mod.lib")
    CROSSLINKS_FILE: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'homology' / "crosslinks.csv")

    def __init__(self, **data):
        super().__init__(**data)
        self.set_mode()

    def set_mode(self):
        if self.sequence_generator and self.geometry_generator:
            self.mode = OperationMode.BOTH
        elif self.sequence_generator:
            self.mode = OperationMode.SEQUENCE
        elif self.geometry_generator:
            self.mode = OperationMode.GEOMETRY
        elif self.fibril:
            self.mode = OperationMode.FIBRIL
        elif self.mix_bool:
            self.mode = OperationMode.MIX
        elif self.replace_bool:
            self.mode = OperationMode.REPLACE

    class Config:
        use_enum_values = True

    def validate_paths(self):
        input_paths = [
            self.PROJECT_ROOT,
            self.DATA_DIR,
            self.HOMOLOGY_LIB_DIR,
            self.TEMPLATE_PDB_PATH,
            self.TEMPLATE_FASTA_PATH,
            self.RESTYP_LIB_PATH,
            self.TOP_HEAV_LIB_PATH,
            self.PAR_MOD_LIB_PATH,
            self.CROSSLINKS_FILE,
            self.file,
            self.crystalcontacts_file,
            self.connect_file,
            self.replace_file,
        ]
        for path in input_paths:
            if isinstance(path, Path) and path is not None and not path.exists():
                raise FileNotFoundError(f"Required input file or directory not found: {path}")

    def update(self, new_config: Dict[str, Any]):
        for key, value in new_config.items():
            if hasattr(self, key):
                setattr(self, key, value)
        self.set_mode()

    def __str__(self):
        return f"ColbuilderConfig(mode={self.mode}, file={self.file}, output={self.output})"

_config_instance = None

def get_config(**kwargs) -> ColbuilderConfig:
    global _config_instance
    if _config_instance is None:
        _config_instance = ColbuilderConfig(**kwargs)
        _config_instance.validate_paths()
    return _config_instance

def load_yaml_config(yaml_file: Path) -> Dict[str, Any]:
    with open(yaml_file, 'r') as file:
        return yaml.safe_load(file)

def validate_config(config: Dict[str, Any]) -> ColbuilderConfig:
    try:
        return ColbuilderConfig(**config)
    except ValidationError as e:
        raise ValueError(f"Configuration validation failed: {e}")