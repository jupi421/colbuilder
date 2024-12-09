"""
Colbuilder Main Geometry Module
Main entry point for geometry operations, coordinating the crystal,
mixing, and replacement services.
"""
from pathlib import Path
from typing import Optional
from ..utils.exceptions import GeometryGenerationError
from ..utils.config import ColbuilderConfig
from ..utils.logger import setup_logger
from .crystal_builder import CrystalBuilder
from .crosslink_mixer import CrosslinkMixer
from .geometry_replacer import GeometryReplacer
from .system import System

LOG = setup_logger(__name__)

class GeometryService:
    """Main service coordinating geometry operations."""
    
    def __init__(self, config: ColbuilderConfig):
        """Initialize geometry service with configuration."""
        self.config = config
        self.crystal_service = CrystalBuilder()
        self.mixer_service = CrosslinkMixer() 
        self.replacer_service = GeometryReplacer()

    async def build_geometry(self) -> Optional[System]:
        """
        Main entry point for geometry generation.
        
        Returns:
            Optional[System]: Generated system or None
        """
        try:
            if not self.config.geometry_generator:
                LOG.info('Set -geometry flag to generate microfibrillar structure PDB file')
                return None

            system = await self.crystal_service.build(self.config)

            if self.config.mix_bool:
                system = await self.mixer_service.mix(system, self.config)

            if self.config.replace_bool:
                system = await self.replacer_service.replace(system, self.config)

            return system

        except GeometryGenerationError:
            raise
        except Exception as e:
            raise GeometryGenerationError(
                message="Unexpected error in geometry generation",
                original_error=e,
                error_code="GEO_ERR_001",
                context={"config": self.config.model_dump()}
            )

async def build_geometry(config: ColbuilderConfig) -> Optional[System]:
    """Build geometry from configuration."""
    service = GeometryService(config)
    return await service.build_geometry()

async def mix_geometry(system: System, config: ColbuilderConfig) -> System:
    """Mix geometry types in system."""
    mixer = CrosslinkMixer()
    return await mixer.mix(system, config)

async def replace_geometry(system: System, config: ColbuilderConfig) -> System:
    """Replace elements in system geometry."""
    replacer = GeometryReplacer()
    return await replacer.replace(system, config)