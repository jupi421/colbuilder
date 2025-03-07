# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from enum import Flag, auto, Enum
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union, Literal, Any
from pydantic import BaseModel, Field, field_validator, computed_field, model_validator, ValidationInfo
import yaml
import re
import os
import subprocess
import json
from functools import lru_cache

from .validators import BioformatValidator
from .exceptions import (
    ConfigurationError,
    SequenceGenerationError,
    GeometryGenerationError,
    TopologyGenerationError
)

from colbuilder.core.utils.logger import setup_logger
LOG = setup_logger(__name__)

class OperationMode(Flag):
    """Operation modes for the Colbuilder pipeline."""
    NONE = 0
    SEQUENCE = auto()
    GEOMETRY = auto()
    TOPOLOGY = auto()
    MIX = auto()
    REPLACE = auto()

class ColbuilderConfig(BaseModel):
    """Main configuration class for the Colbuilder pipeline."""
    
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

    # Fibril geometry generation mode
    geometry_generator: bool = Field(default=False, description="Run geometry generation")
    pdb_file: Optional[Path] = Field(None, description="Input PDB file")
    contact_distance: Optional[float] = Field(None, description="Contact distance for microfibril")
    fibril_length: float = Field(default=334, description="Length of microfibril")
    crystalcontacts_file: Optional[Path] = Field(None, description="Crystal contacts file")
    connect_file: Optional[Path] = Field(None, description="Connect file")
    crystalcontacts_optimize: bool = Field(default=False, description="Optimize crystal contacts")
    solution_space: Union[List[float], Tuple[float, float, float]] = Field(
        default=(1, 1, 1),
        description="Solution space"
    )
    pdb_first_line: Optional[str] = Field(
        default="CRYST1   39.970   26.950  677.900  89.24  94.59 105.58 P 1           2",
        description="Crystal contacts information"
    )
   
    # Mix crosslinks mode
    mix_bool: bool = Field(default=False, description="Generate a mixed crosslinked microfibril")
    ratio_mix: Union[str, Dict[str, int]] = Field(default={}, alias="ratio_mix")
    files_mix: Union[List[Path], Tuple[Path, ...]] = Field(
        default=(),
        description="PDB files with different crosslink types"
    )
    
    # Replace crosslinks mode
    replace_bool: bool = Field(default=False, description="Generate a microfibril with less crosslinks")
    ratio_replace: Optional[float] = Field(None, description="Ratio of crosslinks to be replaced")
    replace_file: Optional[Path] = Field(None, description="File with crosslinks to be replaced")
    
    # Topology generation mode
    topology_generator: bool = Field(default=False, description="Generate topology files")
    force_field: Optional[str] = Field(None, description="Force field to be used")
    martinize2_command: Optional[str] = Field(None, description="Detected Martinize2 command")
    martinize2_env: Optional[str] = Field(None, description="Detected Martinize2 environment")
    use_conda_run: bool = Field(default=False, description="Whether to use conda run for Martinize2")
    go_epsilon: float = Field(default=9.414, description="GO epsilon value for Martini3-CG parametrization (default: 9.414)")

    # Path Configuration
    PROJECT_ROOT: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
    )
    DATA_DIR: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data'
    )
    HOMOLOGY_LIB_DIR: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence'
    )
    TEMPLATE_PDB_PATH: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "template.pdb"
    )
    TEMPLATE_FASTA_PATH: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "template.fasta"
    )
    RESTYP_LIB_PATH: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "modeller" / "restyp_mod.lib"
    )
    TOP_HEAV_LIB_PATH: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "modeller" / "top_heav_mod.lib"
    )
    PAR_MOD_LIB_PATH: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "modeller" / "par_mod.lib"
    )
    CROSSLINKS_FILE: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'sequence' / "crosslinks.csv"
    )
    CHIMERA_SCRIPTS_DIR: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'chimera_scripts'
    )
    FORCE_FIELD_DIR: Path = Field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent / 'data' / 'topology'
    )

    _species_map: Set[str] = {
        "homo_sapiens", "ailuropoda_melanoleuca", "danio_rerio", "mus_musculus",
        "mustela_putorius", "myotis_lucifugus", "otolemur_garnettii",
        "pan_troglodytes", "pongo_abelii", "rattus_norvegicus", "bos_taurus",
        "callithrix_jacchus", "canis_lupus", "loxodonta_africana",
        "oreochromis_niloticus", "oryzias_latipes", "pelodiscus_sinensis",
        "tetraodon_nigroviridis", "xiphophorus_maculatus"
    }

    def __init__(self, **data):
        super().__init__(**data)
        self.ratio_mix = self._convert_ratio_mix(self.ratio_mix)
        self.solution_space = self._convert_to_tuple(self.solution_space)
        self.files_mix = tuple(self.files_mix)
        self.set_mode()

    def model_post_init(self, __context: Any) -> None:
        """Post initialization validation and processing."""
        # Only validate species and handle fasta_file if sequence_generator is True
        if self.sequence_generator:
            if self.fasta_file is None:
                if self.species in self._species_map:
                    fasta_name = f"{self.species.replace('_', '')}.fasta"
                    self.fasta_file = str(self.HOMOLOGY_LIB_DIR / "fasta_sequences" / fasta_name)
                else:
                    error_info = ConfigurationError.get_error_info("CFG_ERR_003")
                    custom_message = f"Must provide fasta_file when using custom species: {self.species}"
                    full_message = f"ERROR: {error_info.message}. {custom_message}"
                    raise ConfigurationError(
                        message=full_message,
                        error_code="CFG_ERR_003"
                    )
        # If sequence_generator is False, validate that PDB file is provided
        elif not self.sequence_generator and self.species not in self._species_map:
            if self.pdb_file is None:
                raise ConfigurationError(
                    "When using a custom species with sequence_generator=False, a PDB file must be provided",
                    error_code="CFG_ERR_003"
                )
        
        self.set_mode()
        
    @model_validator(mode='after')
    def validate_geometry_requirements(self) -> 'ColbuilderConfig':
        """Validate geometry generation requirements."""
        if self.geometry_generator:
            if self.contact_distance is None and self.crystalcontacts_file is None:
                raise ConfigurationError(
                    "Either contact_distance or crystalcontacts_file must be provided for geometry generation",
                    error_code="CFG_ERR_006"
                )
            if self.contact_distance == 0 and self.crystalcontacts_file is None:
                pass
        return self

    @field_validator('contact_distance')
    def validate_contact_distance(cls, value: Optional[float]) -> Optional[float]:
        """Validate contact distance value."""
        if value is not None:
            if not isinstance(value, (int, float)):
                raise ConfigurationError(
                    "Contact distance must be a numeric value",
                    error_code="CFG_ERR_006"
                )
            if value <= 0:
                raise ConfigurationError(
                    f"Contact distance must be positive, got {value}",
                    error_code="CFG_ERR_006"
                )
        return value

    @field_validator('fibril_length')
    def validate_fibril_length(cls, value: float) -> float:
        """Validate fibril length value."""
        if not isinstance(value, (int, float)):
            raise ConfigurationError(
                "Fibril length must be a numeric value",
                error_code="CFG_ERR_006"
            )
        if value <= 0:
            raise ConfigurationError(
                f"Fibril length must be positive, got {value}",
                error_code="CFG_ERR_006"
            )
        if value > 334:
            raise ConfigurationError(
                f"Fibril length must be less than 334, got {value}",
                error_code="CFG_ERR_006"
            )
        return value

    @field_validator('force_field')
    def validate_force_field(cls, value: Optional[str]) -> Optional[str]:
        """Validate force field value."""
        valid_fields = {"amber99", "martini3"}
        if value is not None and value not in valid_fields:
            raise ConfigurationError(
                f"Force field must be one of {valid_fields}, got {value}",
                error_code="CFG_ERR_006"
            )
        return value
    @model_validator(mode='after')
    def validate_species_requirements(self) -> 'ColbuilderConfig':
        """Validate species and related requirements."""
        # If sequence generator is True, species must be from predefined list or fasta_file must be provided
        if self.sequence_generator:
            if self.species not in self._species_map and self.fasta_file is None:
                raise ConfigurationError(
                    f"When sequence_generator is True, species must be one of {self._species_map} or fasta_file must be provided",
                    error_code="CFG_ERR_003"
                )
        # If sequence generator is False and species is not from predefined list, PDB file must be provided
        elif self.species not in self._species_map:
            if self.pdb_file is None:
                raise ConfigurationError(
                    "When using a custom species with sequence_generator=False, a PDB file must be provided",
                    error_code="CFG_ERR_003"
                )
        return self

    @field_validator('species', mode='before')
    def convert_species_to_lowercase(cls, value: str) -> str:
        """Convert species name to lowercase."""
        return value.lower() if value else value

    def _convert_ratio_mix(self, value: Union[str, Dict[str, int]]) -> Dict[str, int]:
        """Convert string ratio mix to dictionary."""
        if isinstance(value, str):
            return {item.split(':')[0]: int(item.split(':')[1]) for item in value.split()}
        return value

    def _convert_to_tuple(self, value: Union[List[float], Tuple[float, float, float]]) -> Tuple[float, float, float]:
        """Convert list to tuple for solution space."""
        if isinstance(value, list):
            return tuple(value)
        return value

    @field_validator('solution_space', mode='before')
    def validate_solution_space(cls, value):
        """Validate solution space format."""
        if isinstance(value, list):
            return tuple(value)
        return value

    @field_validator('files_mix', mode='before')
    def validate_files_mix(cls, value):
        """Validate files mix format."""
        return tuple(value)

    @property
    def ratio_mix(self) -> Dict[str, int]:
        """Get ratio mix property."""
        return self._ratio_mix

    @ratio_mix.setter
    def ratio_mix(self, value: Union[str, Dict[str, int]]):
        """Set ratio mix property."""
        self._ratio_mix = self._convert_ratio_mix(value)

    def set_mode(self):
        """Set operation mode based on configuration flags."""
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

    def validate_paths(self):
        """Validate existence of required input paths and files."""
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
                raise ConfigurationError(
                    f"Required input file or directory not found: {path}",
                    error_code="CFG_ERR_002"
                )

    @field_validator('ratio_mix', mode='before')
    def validate_ratio_mix(cls, v: Union[str, Dict[str, int]], info: ValidationInfo) -> Dict[str, int]:
        """Validate ratio mix format and values."""
        try:
            if isinstance(v, str):
                result = {item.split(':')[0]: int(item.split(':')[1]) for item in v.split()}
                if sum(result.values()) != 100:
                    raise ConfigurationError(
                        "Mix ratios must sum to 100%",
                        error_code="CFG_ERR_004"
                    )
                return result
            return v
        except (ValueError, IndexError):
            raise ConfigurationError(
                "Invalid ratio mix format. Expected 'Type:percentage Type:percentage'",
                error_code="CFG_ERR_004"
            )

    @model_validator(mode='after')
    def validate_mix_config(self) -> 'ColbuilderConfig':
        """Validate mixing configuration."""
        if self.mix_bool:
            if not self.ratio_mix:
                raise ConfigurationError(
                    "ratio_mix must be specified when mix_bool is True",
                    error_code="CFG_ERR_004"
                )
            if not self.files_mix:
                raise ConfigurationError(
                    "files_mix must be specified when mix_bool is True",
                    error_code="CFG_ERR_004"
                )

        if self.n_term_combination:
            if not re.match(r'^\d+\.[A-C]\s*-\s*\d+\.[A-C]$', self.n_term_combination):
                raise ConfigurationError(
                    f"Invalid N-terminal combination format: {self.n_term_combination}. "
                    "Expected format: 'ResidueNumber.Chain - ResidueNumber.Chain' (e.g., '9.C - 947.A')",
                    error_code="CFG_ERR_006"
                )
    
        if self.c_term_combination:
            if not re.match(r'^\d+\.[A-C]\s*-\s*\d+\.[A-C]$', self.c_term_combination):
                raise ConfigurationError(
                    f"Invalid C-terminal combination format: {self.c_term_combination}. "
                    "Expected format: 'number.Chain - number.Chain' (e.g., '1047.C - 104.C')",
                    error_code="CFG_ERR_006"
                )

        if self.replace_bool:
            if self.ratio_replace is None:
                raise ConfigurationError(
                    "ratio_replace must be specified when replace_bool is True",
                    error_code="CFG_ERR_006"
                )
            if not isinstance(self.ratio_replace, (int, float)):
                raise ConfigurationError(
                    "ratio_replace must be a numeric value",
                    error_code="CFG_ERR_006"
                )
            if self.ratio_replace < 0 or self.ratio_replace > 100:
                raise ConfigurationError(
                    f"ratio_replace must be between 0 and 100, got {self.ratio_replace}",
                    error_code="CFG_ERR_006"
                )

        if len(self.solution_space) != 3:
            raise ConfigurationError(
                f"Solution space must have exactly 3 values, got {len(self.solution_space)}",
                error_code="CFG_ERR_006"
            )
        for val in self.solution_space:
            if not isinstance(val, (int, float)) or val <= 0:
                raise ConfigurationError(
                    f"Solution space values must be positive numbers, got {val}",
                    error_code="CFG_ERR_006"
                )

        return self
    
    def get_conda_env_path(self) -> Optional[str]:
        """
        Get the full path to the conda environment for martinize2.
        
        Returns:
            Optional[str]: Full path to the conda environment or None if not found
        """
        if not self.martinize2_env:
            return None
            
        try:
            # Run conda info command to get environment paths
            result = subprocess.run(
                ["conda", "info", "--envs", "--json"],
                capture_output=True, text=True, check=True
            )
            env_data = json.loads(result.stdout)
            
            # Look for the specified environment in the paths
            for env_path in env_data["envs"]:
                if os.path.basename(env_path) == self.martinize2_env:
                    LOG.debug(f"Found conda environment path: {env_path}")
                    return env_path
            
            # Try alternative method using CONDA_PREFIX
            if "CONDA_PREFIX" in os.environ:
                base_path = os.environ["CONDA_PREFIX"]
                # Go up one level if we're in an environment already
                if os.path.basename(base_path) not in ["anaconda3", "miniconda3"]:
                    base_path = os.path.dirname(base_path)
                candidate_path = os.path.join(os.path.dirname(base_path), "envs", self.martinize2_env)
                if os.path.exists(candidate_path):
                    LOG.debug(f"Found conda environment path via CONDA_PREFIX: {candidate_path}")
                    return candidate_path
            
            LOG.warning(f"Could not find conda environment path for: {self.martinize2_env}")
            return None
        except Exception as e:
            LOG.error(f"Error determining conda environment path: {e}")
            return None

    def update(self, new_config: Dict[str, Any]):
        """Update configuration with new values."""
        for key, value in new_config.items():
            if hasattr(self, key):
                setattr(self, key, value)
        self.set_mode()

    def __str__(self):
        """String representation of the configuration."""
        return f"ColbuilderConfig(mode={self.mode}, pdb_file={self.pdb_file}, output={self.output})"
    
    @property
    def output(self) -> str:
        """Get the output filename."""
        return f"collagen_fibril_{self.species}"

    class Config:
        """Pydantic configuration."""
        use_enum_values = True


def validate_input_files(config: ColbuilderConfig) -> None:
    """Validate input files if they are provided and required."""
    validator = BioformatValidator()
    
    try:
        if config.sequence_generator and config.fasta_file is not None:
            fasta_path = Path(config.fasta_file)
            warnings = validator.validate_input_files(fasta_path, None)
            if warnings:
                for warning in warnings:
                    LOG.warning(warning)
                    
        if config.geometry_generator and config.pdb_file is not None:
            pdb_path = config.pdb_file
            warnings = validator.validate_input_files(
                fasta_path=Path(config.fasta_file) if config.fasta_file else None, 
                pdb_path=pdb_path
            )
            if warnings:
                for warning in warnings:
                    LOG.warning(warning)
                        
    except SequenceGenerationError as e:
        raise ConfigurationError(
            f"ERROR: {str(e)}",
            original_error=e,
            error_code="CFG_ERR_005"
        )
    except GeometryGenerationError as e:
        raise ConfigurationError(
            f"ERROR: {str(e)}",
            original_error=e,
            error_code="CFG_ERR_005"
        )
    except Exception as e:
        raise ConfigurationError(
            "File validation failed",
            original_error=e,
            error_code="CFG_ERR_005"
        )


_config_instance = None

def get_config(**kwargs) -> ColbuilderConfig:
    """Get or create the configuration singleton."""
    global _config_instance
    if _config_instance is None:
        try:
            _config_instance = ColbuilderConfig(**kwargs)
            _config_instance.validate_paths()
            validate_input_files(_config_instance)
            
            # Detect Martinize2 if topology generation is enabled
            if _config_instance.topology_generator:
                try:
                    from colbuilder.core.utils.martinize_finder import find_martinize2_executable
                    executable_path, use_conda, conda_env = find_martinize2_executable()
                    _config_instance.martinize2_command = executable_path
                    _config_instance.martinize2_env = conda_env
                    _config_instance.use_conda_run = use_conda
                    LOG.info(f"Detected Martinize2: {executable_path}" + 
                            (f" (using conda environment: {conda_env})" if use_conda else ""))
                except Exception as e:
                    LOG.warning(f"Failed to detect Martinize2: {str(e)}")
                    LOG.warning("Topology generation may fail if Martinize2 is not available")
                    
        except ConfigurationError:
            raise
        except Exception as e:
            raise ConfigurationError(
                "Configuration initialization failed",
                original_error=e,
                error_code="CFG_ERR_001",
                context={"kwargs": kwargs} 
            )
    return _config_instance


def load_yaml_config(yaml_file: Path) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    try:
        with open(yaml_file, 'r') as file:
            return yaml.safe_load(file)
    except Exception as e:
        raise ConfigurationError(
            f"ERROR: {yaml_file}",
            original_error=e,
            error_code="CFG_ERR_001"
        )


def validate_config(config: Dict[str, Any]) -> ColbuilderConfig:
    """Validate configuration dictionary and create config object."""
    try:
        return ColbuilderConfig(**config)
    except Exception as e:
        raise ConfigurationError(
            "Configuration validation failed",
            original_error=e,
            error_code="CFG_ERR_001"
        )
