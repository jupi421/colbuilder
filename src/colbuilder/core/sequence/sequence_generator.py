# sequence_generator.py
import asyncio
import os
import tempfile
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple, List, Optional, Dict, Any
from contextlib import asynccontextmanager
from dataclasses import asdict

from colbuilder.core.utils.constants import (
    TEMP_FILE_SUFFIX,
    DISORIENTED_SUFFIX,
    PDB_EXTENSION,
    FASTA_EXTENSION
)
from colbuilder.core.utils.data_structures import CrosslinkPair, OptimizationState
from colbuilder.core.utils.crosslinks import (
    CrosslinkOptimizer,
    extract_crosslinks_from_dataframe
)
from colbuilder.core.utils.files import (
    update_pdb_header,
    suppress_output
)

from colbuilder.core.sequence.alignment import align_sequences
from colbuilder.core.sequence.modeller import run_modeller
from colbuilder.core.sequence.mutate_crosslinks import apply_crosslinks
from colbuilder.core.sequence.optimize_crosslinks import optimize_structure

from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.exceptions import SequenceGenerationError, SystemError
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

class SequenceGenerator:
    """
    Manages the generation of collagen structure from sequences.
    
    This class coordinates the entire structure from sequence generation process, including:
    - Sequence alignment
    - Structure modeling
    - Crosslink application and optimization
    - File management and cleanup
    
    Attributes:
        config (ColbuilderConfig): Configuration for sequence generation
        _temp_dir (Optional[Path]): Path to temporary working directory
        _crosslinks (List[CrosslinkPair]): List of crosslinks to apply
        _state (Dict[str, Any]): Current generation state
    """
    
    def __init__(self, config: ColbuilderConfig):
        """
        Initialize the structure from sequence generator.
        
        Args:
            config: Configuration object for structure generation
        """
        self.config = config
        self._temp_dir: Optional[Path] = None
        self._crosslinks: List[CrosslinkPair] = []
        self._state: Dict[str, Any] = {}
        self.steps = 4 if config.crosslink else 2
        
    @asynccontextmanager
    async def _manage_resources(self):
        """Context manager for resource lifecycle management."""
        try:
            original_dir = Path.cwd()
            working_dir = Path(self.config.working_directory)

            if not working_dir.exists():
                LOG.error(f"Working directory doesn't exist: {working_dir}")
                raise SequenceGenerationError(
                    "Working directory does not exist",
                    error_code="SEQ_ERR_004",
                    context={"working_directory": str(working_dir)}
                )

            temp_dir = (
                working_dir / "temp_working_dir"
                if self.config.debug
                else Path(tempfile.mkdtemp(dir=working_dir))
            )
            temp_dir = temp_dir.resolve()  
            temp_dir.mkdir(exist_ok=True)
            self._temp_dir = temp_dir
            os.chdir(temp_dir)

            yield

        except Exception as e:
            LOG.error(f"Error setting up directories: {e}")
            raise SequenceGenerationError(
                "Failed to set up working environment",
                error_code="SEQ_ERR_004",
                context={
                    "working_dir": str(working_dir),
                    "temp_dir": str(temp_dir) if 'temp_dir' in locals() else None,
                    "current_dir": str(Path.cwd()),
                    "error": str(e),
                },
            )
        finally:
            try:
                os.chdir(original_dir)
            except Exception as e:
                LOG.warning(f"Failed to restore original directory {original_dir}: {e}")

            if not self.config.debug and self._temp_dir and self._temp_dir.exists():
                try:
                    shutil.rmtree(self._temp_dir)
                except Exception as e:
                    LOG.warning(f"Failed to clean up temporary directory {self._temp_dir}: {e}")

            self._temp_dir = None

                
    async def _load_crosslinks(self) -> None:
        """
        Load crosslink information from configuration.
            
        Raises:
            SequenceGenerationError: If crosslink loading fails
        """
        try:
            crosslinks_df = pd.read_csv(self.config.CROSSLINKS_FILE)
            species_crosslinks = crosslinks_df[crosslinks_df['species'] == self.config.species]
                
            if species_crosslinks.empty:
                raise SequenceGenerationError(
                    f"No crosslinks found for species: {self.config.species}",
                    error_code="SEQ_ERR_002",
                    context={"species": self.config.species}
                )
                    
            self._crosslinks = []
                
            n_crosslinks = extract_crosslinks_from_dataframe(
                species_crosslinks,
                "N",
                self.config.n_term_type,
                self.config.n_term_combination
            )
            self._crosslinks.extend(n_crosslinks)
                
            c_crosslinks = extract_crosslinks_from_dataframe(
                species_crosslinks,
                "C",
                self.config.c_term_type,
                self.config.c_term_combination
            )
            self._crosslinks.extend(c_crosslinks)
                
        except SequenceGenerationError:
            raise
        except Exception as e:
            raise SequenceGenerationError(
                "Error loading crosslinks",
                original_error=e,
                error_code="SEQ_ERR_002",
                context={
                    "crosslinks_file": str(self.config.CROSSLINKS_FILE),
                    "species": self.config.species
                }
            )
            
    async def _run_alignment(self, fasta_path: Path) -> Tuple[Path, Path]:
        """Run sequence alignment."""
        try:
            LOG.info(f"Step 1/{self.steps} Performing sequence alignment")
            
            fasta_path = fasta_path.resolve()
            template_fasta = Path(self.config.TEMPLATE_FASTA_PATH).resolve()
            template_pdb = Path(self.config.TEMPLATE_PDB_PATH).resolve()
            
            file_prefix = fasta_path.stem
            LOG.debug(f"Template FASTA: {template_fasta}")
            LOG.debug(f"Template PDB: {template_pdb}")
            
            if not fasta_path.exists():
                raise SequenceGenerationError(
                    "Input FASTA file not found",
                    error_code="SEQ_ERR_004",
                    context={"fasta_path": str(fasta_path)}
                )
                
            if not template_fasta.exists():
                raise SequenceGenerationError(
                    "Template FASTA file not found",
                    error_code="SEQ_ERR_004",
                    context={"template_fasta": str(template_fasta)}
                )
                
            if not template_pdb.exists():
                raise SequenceGenerationError(
                    "Template PDB file not found",
                    error_code="SEQ_ERR_004",
                    context={"template_pdb": str(template_pdb)}
                )
            
            with suppress_output():
                msa_output_path, modeller_output = align_sequences(
                    fasta_path,
                    template_fasta,
                    file_prefix,
                    template_pdb
                )
                
            msa_output = Path(msa_output_path).resolve()
            modeller_out = Path(modeller_output).resolve()
            
            LOG.debug(f"Alignment completed. MSA output: {msa_output}")
            LOG.debug(f"Modeller output: {modeller_out}")
            
            if not msa_output.exists():
                raise SequenceGenerationError(
                    "MSA output file not created",
                    error_code="SEQ_ERR_001",
                    context={"msa_output": str(msa_output)}
                )
                
            if not modeller_out.exists():
                raise SequenceGenerationError(
                    "Modeller output file not created",
                    error_code="SEQ_ERR_001",
                    context={"modeller_output": str(modeller_out)}
                )
                
            return msa_output, modeller_out
            
        except Exception as e:
            LOG.error(f"Sequence alignment failed: {str(e)}", exc_info=True)
            raise SequenceGenerationError(
                "Sequence alignment failed",
                error_code="SEQ_ERR_001",
                context={
                    "fasta_file": str(fasta_path),
                    "template_fasta": str(self.config.TEMPLATE_FASTA_PATH),
                    "current_dir": str(Path.cwd()),
                    "error": str(e)
                }
            )

    async def _run_modelling(self, modeller_output: Path, file_prefix: str) -> Path:
        """Run MODELLER for structure generation."""
        try:
            LOG.info(f"Step 2/{self.steps} Generating structure with MODELLER")
            
            modeller_output = modeller_output.resolve()
            template_pdb = Path(self.config.TEMPLATE_PDB_PATH).resolve()
            restyp_lib = Path(self.config.RESTYP_LIB_PATH).resolve()
            top_heav_lib = Path(self.config.TOP_HEAV_LIB_PATH).resolve()
            par_mod_lib = Path(self.config.PAR_MOD_LIB_PATH).resolve()
            
            LOG.debug(f"Using libraries:")
            LOG.debug(f"  Template PDB: {template_pdb}")
            LOG.debug(f"  RESTYP: {restyp_lib}")
            LOG.debug(f"  TOP_HEAV: {top_heav_lib}")
            LOG.debug(f"  PAR_MOD: {par_mod_lib}")
            
            for path, name in [
                (modeller_output, "Modeller input"),
                (template_pdb, "Template PDB"),
                (restyp_lib, "RESTYP library"),
                (top_heav_lib, "TOP_HEAV library"),
                (par_mod_lib, "PAR_MOD library")
            ]:
                if not path.exists():
                    raise SequenceGenerationError(
                        f"{name} file not found",
                        error_code="SEQ_ERR_004",
                        context={"missing_file": str(path)}
                    )
            
            with suppress_output(): 
                output_pdb = run_modeller(
                    aligned_file=str(modeller_output),
                    template_pdb=str(template_pdb),
                    output_prefix=file_prefix,
                    restyp_lib=str(restyp_lib),
                    top_heav_lib=str(top_heav_lib),
                    par_mod_lib=str(par_mod_lib)
                )
                
            output_path = Path(output_pdb).resolve()
            LOG.debug(f"MODELLER completed. Output PDB: {output_path}")
            
            if not output_path.exists():
                raise SequenceGenerationError(
                    "MODELLER output file not created",
                    error_code="SEQ_ERR_001",
                    context={"output_pdb": str(output_path)}
                )
                
            return output_path
            
        except Exception as e:
            LOG.error(f"MODELLER failed: {str(e)}", exc_info=True)
            raise SequenceGenerationError(
                "MODELLER structure generation failed",
                error_code="SEQ_ERR_001",
                context={
                    "modeller_output": str(modeller_output),
                    "template_pdb": str(self.config.TEMPLATE_PDB_PATH),
                    "current_dir": str(Path.cwd()),
                    "error": str(e)
                }
            )

    async def generate(self) -> Tuple[Path, Path]:
        """Generate collagen structure from sequence."""
        LOG.debug(f"Working directory: {self.config.working_directory}")
        LOG.debug(f"Input FASTA: {self.config.fasta_file}")
        LOG.debug(f"Debug mode: {self.config.debug}")

        original_dir = Path.cwd().resolve()
        
        async with self._manage_resources():
            try:
                fasta_file = Path(self.config.fasta_file).resolve()
                
                if not fasta_file.exists():
                    LOG.error(f"FASTA file not found: {fasta_file}")
                    raise SequenceGenerationError(
                        "FASTA file not found",
                        error_code="SEQ_ERR_004",
                        context={"file_path": str(fasta_file)}
                    )
                
                file_prefix = fasta_file.stem
                fasta_copy = self._temp_dir / fasta_file.name
                fasta_copy.write_bytes(fasta_file.read_bytes())
                
                msa_output, modeller_output = await self._run_alignment(fasta_copy)
                LOG.debug(f"Alignment complete - MSA: {msa_output}, Modeller: {modeller_output}")
                
                output_pdb = await self._run_modelling(modeller_output, file_prefix)
                LOG.debug(f"MODELLER complete - Output: {output_pdb}")
                
                if self.config.crosslink:
                    await self._load_crosslinks()
                    with suppress_output():
                        output_pdb = await self._apply_crosslinks(output_pdb, file_prefix)
                        LOG.debug(f"Crosslinks applied - Output: {output_pdb}")
                
                update_pdb_header(output_pdb, str(self.config.pdb_first_line))
                
                final_output = await self._finalize_output(output_pdb, file_prefix)
                
                final_msa = original_dir / msa_output.name
                final_pdb = original_dir / final_output.name

                msa_output = msa_output.resolve()
                final_output = final_output.resolve()
                
                try:
                    shutil.copy2(msa_output, final_msa)
                    shutil.copy2(final_output, final_pdb)
                    
                    if not final_msa.exists():
                        raise SequenceGenerationError(
                            "Failed to copy MSA file to working directory",
                            error_code="SEQ_ERR_004",
                            context={"source": str(msa_output), "destination": str(final_msa)}
                        )
                    if not final_pdb.exists():
                        raise SequenceGenerationError(
                            "Failed to copy PDB file to working directory",
                            error_code="SEQ_ERR_004",
                            context={"source": str(final_output), "destination": str(final_pdb)}
                        )
                        
                except Exception as e:
                    LOG.error(f"Error copying output files: {e}")
                    raise SequenceGenerationError(
                        "Failed to copy output files",
                        error_code="SEQ_ERR_004",
                        context={
                            "msa_source": str(msa_output),
                            "pdb_source": str(final_output),
                            "working_dir": str(original_dir),
                            "error": str(e)
                        }
                    )

                return final_msa, final_pdb
                
            except (SequenceGenerationError, SystemError):
                LOG.error("Known error occurred during generation", exc_info=True)
                raise
            except Exception as e:
                LOG.error(f"Unexpected error during generation: {str(e)}", exc_info=True)
                raise SequenceGenerationError(
                    "Unexpected error during sequence generation",
                    error_code="SEQ_ERR_001",
                    context={
                        "state": self._state,
                        "error": str(e),
                        "error_type": type(e).__name__
                    }
                )

    async def _apply_crosslinks(self, input_pdb: Path, file_prefix: str) -> Path:
        """
        Apply crosslinks to the PDB structure.
        
        Args:
            input_pdb: Path to input PDB file
            file_prefix: Prefix for output files
            
        Returns:
            Path to PDB file with crosslinks applied
            
        Raises:
            SequenceGenerationError: If crosslink application fails
        """
        try:
            n_suffix = f"N_{self.config.n_term_type}" if self.config.n_term_type else "N_NONE"
            c_suffix = f"C_{self.config.c_term_type}" if self.config.c_term_type else "C_NONE"
            output_pdb_crosslinked = f"{file_prefix}_{n_suffix}_{c_suffix}_temp.pdb"
            
            LOG.info(f"Step 3/{self.steps} Applying crosslinks to structure")
            
            n_crosslinks = [c for c in self._crosslinks if c.terminal_type == 'N']
            c_crosslinks = [c for c in self._crosslinks if c.terminal_type == 'C']
            
            n_crosslink_data = None if not n_crosslinks else {
                'P1': n_crosslinks[0].position1.position_str,
                'R1': n_crosslinks[0].position1.residue_type,
                'A1': n_crosslinks[0].position1.atom_name,
                'P2': n_crosslinks[0].position2.position_str,
                'R2': n_crosslinks[0].position2.residue_type,
                'A2': n_crosslinks[0].position2.atom_name
            }
            
            c_crosslink_data = None if not c_crosslinks else {
                'P1': c_crosslinks[0].position1.position_str,
                'R1': c_crosslinks[0].position1.residue_type,
                'A1': c_crosslinks[0].position1.atom_name,
                'P2': c_crosslinks[0].position2.position_str,
                'R2': c_crosslinks[0].position2.residue_type,
                'A2': c_crosslinks[0].position2.atom_name
            }
            
            result_path = Path(apply_crosslinks(
                str(input_pdb),
                output_pdb_crosslinked,
                n_crosslink_data,
                c_crosslink_data,
                self.config
            ))
            
            if not result_path.exists():
                raise SequenceGenerationError(
                    "Crosslink application failed to create output file",
                    error_code="SEQ_ERR_002",
                    context={"output_path": str(result_path)}
                )
                
            return result_path
            
        except SequenceGenerationError:
            raise
        except Exception as e:
            raise SequenceGenerationError(
                "Error applying crosslinks",
                original_error=e,
                error_code="SEQ_ERR_002",
                context={
                    "input_pdb": str(input_pdb),
                    "file_prefix": file_prefix
                }
            )
            
    async def _finalize_output(self, input_pdb: Path, file_prefix: str) -> Path:
        """Finalize output files and perform crosslink optimization if needed."""
        try:
            dis_output_name = input_pdb.stem.replace('_temp', '_disoriented') + input_pdb.suffix
            dis_output_path = self._temp_dir / dis_output_name
            dis_output_path.write_bytes(input_pdb.read_bytes())
            
            final_output_name = input_pdb.stem.replace('_temp', '') + input_pdb.suffix
            temp_final_output = self._temp_dir / final_output_name  # Temporary location
            working_dir_output = Path(self.config.working_directory) / final_output_name  # Final location
            
            if self.config.crosslink and self._crosslinks:
                
                LOG.info(f"Step 4/{self.steps} Optimizing crosslink positions")
                optimizer = CrosslinkOptimizer(
                    crosslink_pairs=self._crosslinks,
                    chimera_scripts_dir=Path(self.config.CHIMERA_SCRIPTS_DIR)
                )
                
                final_distance, final_output = await optimizer.optimize(
                    input_pdb=dis_output_path,
                    output_pdb=temp_final_output
                )
                
                working_dir_output.write_bytes(final_output.read_bytes())
                
                self._state.update({
                    'optimization_distance': final_distance,
                    'optimization_attempts': optimizer.state.attempt_number
                })
                
            else:
                LOG.info("No optimization needed, preparing final output")
                temp_final_output.write_bytes(dis_output_path.read_bytes())
                working_dir_output.write_bytes(temp_final_output.read_bytes())
            
            update_pdb_header(working_dir_output, str(self.config.pdb_first_line))
            
            try:
                dis_output_path.unlink()
            except Exception as e:
                LOG.warning(f"Failed to delete temporary file {dis_output_path}: {str(e)}")
            
            if not working_dir_output.exists():
                raise SequenceGenerationError(
                    "Failed to create final output file",
                    error_code="SEQ_ERR_004",
                    context={"final_output": str(working_dir_output)}
                )
            
            return working_dir_output
            
        except SequenceGenerationError:
            raise
        except Exception as e:
            LOG.error(f"Error in finalize_output: {str(e)}", exc_info=True)
            raise SequenceGenerationError(
                "Error finalizing output",
                error_code="SEQ_ERR_001",
                context={
                    "input_pdb": str(input_pdb),
                    "file_prefix": file_prefix,
                    "state": self._state
                }
            )
            
    @property
    def state(self) -> Dict[str, Any]:
        """Get current generation state."""
        return self._state.copy()
        
    @property
    def has_crosslinks(self) -> bool:
        """Check if structure has crosslinks."""
        return bool(self._crosslinks)
        
    def get_crosslink_info(self) -> List[Dict[str, Any]]:
        """Get information about current crosslinks."""
        return [asdict(crosslink) for crosslink in self._crosslinks]