# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from typing import Any, Dict, Optional, Tuple, List, Set, Union
from pathlib import Path
from pydantic import BaseModel, Field, field_validator, validator
from enum import Flag, auto, Enum
import yaml

class OperationMode(Flag):
    NONE = 0
    SEQUENCE = auto()
    GEOMETRY = auto()
    TOPOLOGY = auto()
    MIX = auto()
    REPLACE = auto()
    
class ColbuilderConfig(BaseModel):
    # Operation mode 
    mode: Optional[OperationMode] = Field(None, description="Operation mode")
    debug: bool = Field(default=False, description="Enable debug logging")
    working_directory: Path = Field(default=Path.cwd(), description="Working directory")
    config_file: Optional[Path] = Field(None, description="YAML configuration file")
    
    # PDB file generation mode (from sequence)
    sequence_generator: bool = Field(default=False, description="Run sequence generation")
    species: str = Field(..., description="Species name")
    fasta_file: Optional[str] = Field(None, description="FASTA file")
    crosslink: bool = Field(default=False, description="Apply crosslinks")
    n_term_type: Optional[str] = Field(None, description="N-terminal type")
    c_term_type: Optional[str] = Field(None, description="C-terminal type")
    n_term_combination: Optional[str] = Field(None, description="N-terminal combination")
    c_term_combination: Optional[str] = Field(None, description="C-terminal combination")

    @field_validator('species', mode='before')
    def convert_species_to_lowercase(cls, value: str) -> str:
        return value.lower() if value else value

    # Fibril geometry generation mode (from pdb file)
    geometry_generator: bool = Field(default=False, description="Run geometry generation")
    pdb_file: Optional[Path] = Field(None, description="Input PDB file")
    contact_distance: Optional[float] = Field(None, description="Contact distance for microfibril")
    fibril_length: float = Field(default=334, description="Length of microfibril")
    crystalcontacts_file: Optional[Path] = Field(None, description="Crystal contacts file")
    connect_file: Optional[Path] = Field(None, description="Connect file")
    crystalcontacts_optimize: bool = Field(default=False, description="Optimize crystal contacts")
    solution_space: Union[List[float], Tuple[float, float, float]] = Field(default=(1, 1, 1), description="Solution space")
    pdb_first_line: Optional[str] = Field(default="CRYST1   39.970   26.950  677.900  89.24  94.59 105.58 P 1           2", description="Crystal contacts information")
    
    # Mix crosslinks mode
    mix_bool: bool = Field(default=False, description="Generate a mixed crosslinked microfibril")
    ratio_mix: Union[str, Dict[str, int]] = Field(default={}, alias="ratio_mix")
    files_mix: Union[List[Path], Tuple[Path, ...]] = Field(default=(), description="PDB files with different crosslink types")
    
    # Replace crosslinks by original residue mode
    replace_bool: bool = Field(default=False, description="Generate a microfibril with less crosslinks")
    ratio_replace: Optional[float] = Field(None, description="Ratio of crosslinks to be replaced")
    replace_file: Optional[Path] = Field(None, description="File with crosslinks to be replaced")
    
    # Topology generation mode
    topology_generator: bool = Field(default=False, description="Generate topology files")
    force_field: Optional[str] = Field(None, description="Force field to be used")
    
    # Path Configuration
    PROJECT_ROOT: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent)
    DATA_DIR: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data')
    HOMOLOGY_LIB_DIR: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence')
    TEMPLATE_PDB_PATH: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "template.pdb")
    TEMPLATE_FASTA_PATH: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "template.fasta")
    RESTYP_LIB_PATH: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "modeller" / "restyp_mod.lib")
    TOP_HEAV_LIB_PATH: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "modeller" / "top_heav_mod.lib")
    PAR_MOD_LIB_PATH: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "modeller" / "par_mod.lib")
    CROSSLINKS_FILE: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "crosslinks.csv")
    CHIMERA_SCRIPTS_DIR: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'chimera_scripts', description="Directory with scripts run as subprocesses by Chimera")
    FORCE_FIELD_DIR: Path = Field(default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'topology', description="Directory containing force field files")

    def __init__(self, **data):
        super().__init__(**data)
        self.ratio_mix = self._convert_ratio_mix(self.ratio_mix)
        self.solution_space = self._convert_to_tuple(self.solution_space)
        self.files_mix = tuple(self.files_mix)
        self.set_mode()

    # Species mapping
    _species_map = {
        "homo_sapiens", "ailuropoda_melanoleuca", "danio_rerio", "mus_musculus",
        "mustela_putorius", "myotis_lucifugus", "otolemur_garnettii",
        "pan_troglodytes", "pongo_abelii", "rattus_norvegicus", "bos_taurus",
        "callithrix_jacchus", "canis_lupus", "loxodonta_africana",
        "oreochromis_niloticus", "oryzias_latipes", "pelodiscus_sinensis",
        "tetraodon_nigroviridis", "xiphophorus_maculatus"
    }
    
    def model_post_init(self, __context: Any) -> None:
        if self.fasta_file is None:
            if self.species in self._species_map:
                fasta_name = f"{self.species.replace('_', '')}.fasta"
                self.fasta_file = str(self.HOMOLOGY_LIB_DIR / "fasta_sequences" / fasta_name)
            else:
                raise ValueError(f"Must provide fasta_file when using custom species: {self.species}")
        self.set_mode()
        
    def _convert_ratio_mix(self, value: Union[str, Dict[str, int]]) -> Dict[str, int]:
        if isinstance(value, str):
            return {item.split(':')[0]: int(item.split(':')[1]) for item in value.split()}
        return value

    def _convert_to_tuple(self, value: Union[List[float], Tuple[float, float, float]]) -> Tuple[float, float, float]:
        if isinstance(value, list):
            return tuple(value)
        return value

    @field_validator('solution_space', mode='before')
    def validate_solution_space(cls, value):
        if isinstance(value, list):
            return tuple(value)
        return value

    @field_validator('files_mix', mode='before')
    def validate_files_mix(cls, value):
        return tuple(value)

    @property
    def ratio_mix(self) -> Dict[str, int]:
        return self._ratio_mix

    @ratio_mix.setter
    def ratio_mix(self, value: Union[str, Dict[str, int]]):
        self._ratio_mix = self._convert_ratio_mix(value)

    def set_mode(self):
        self.mode = OperationMode.NONE
        if self.sequence_generator:
            self.mode |= OperationMode.SEQUENCE
        if self.geometry_generator:
            self.mode |= OperationMode.GEOMETRY
        if self.topology_generator:
            self.mode |= OperationMode.TOPOLOGY
        if self.mix_bool:
            self.mode |= OperationMode.MIX
        if self.replace_bool:
            self.mode |= OperationMode.REPLACE

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
            self.FORCE_FIELD_DIR,
            self.pdb_file,
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
        return f"ColbuilderConfig(mode={self.mode}, pdb_file={self.pdb_file}, output={self.output})"
    
    @property
    def output(self) -> str:
        return f"collagen_fibril_{self.species}"
    
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