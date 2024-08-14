# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import numpy as np
from typing import Optional, Dict, List, Any

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

class Replace:
    """
    Remove crosslinks from certain models to reduce the overall number of crosslinks
    """
    def __init__(self, ratio_replace: Optional[float] = None, system: Optional[Any] = None, fibril_length: Optional[float] = None):
        self.ratio = ratio_replace
        self.system = system
        self.is_line: tuple = ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')
        self.external: Dict = {}
        if system and fibril_length:
            cog = self.system.get_model(model_id=0.0).get_cog()
            if isinstance(cog, np.ndarray):
                cog = cog.tolist()
            self.z_min = cog[2] - fibril_length/2
            self.z_max = cog[2] + fibril_length/2
        LOG.debug(f"Replace initialized with ratio: {ratio_replace}, fibril length: {fibril_length}")
        LOG.debug(f"Z range: {self.z_min} to {self.z_max}")

    def run_replace(self, ratio_replace: Optional[float] = None, system: Optional[Any] = None) -> Any:
        """
        Select models with crosslinks to be replaced by lysines
        Args:
            ratio_replace: Ratio of crosslinks to replace
            system: System object to modify
        Returns:
            Modified system
        """
        if system is None:
            system = self.system
        if ratio_replace is not None:
            self.ratio = ratio_replace
        
        target_ratio = int(self.ratio) / 100
        
        LOG.debug(f"Starting replace process with target ratio: {target_ratio}")
        
        iteration = 0
        max_iterations = 1000000
        while True:
            replace_count = system.count_states(state='replace')
            total_count = system.count_states(state='none') + system.count_states(state='protect') + replace_count
            current_ratio = replace_count / total_count if total_count > 0 else 0
            
            LOG.debug(f"Iteration {iteration}: Current replace ratio: {current_ratio:.4f}, Target ratio: {target_ratio:.4f}")
            LOG.debug(f"Replace count: {replace_count}, Total count: {total_count}")
            
            if current_ratio >= target_ratio:
                break
            
            model_id = self.draw_model(system=system)
            LOG.debug(f"Selected model ID: {model_id}")
            model = system.get_model(model_id=model_id)
            self.replace_model(model=model)
            
            iteration += 1
            if iteration > max_iterations:
                LOG.warning(f"Exceeded maximum iterations ({max_iterations}). Breaking loop.")
                LOG.warning(f"Final replace ratio: {current_ratio:.4f}")
                break
        
        LOG.debug(f"Final replace ratio: {current_ratio:.4f}")
        return system, current_ratio

    def replace_model(self, model: Any) -> None:
        """
        Change status of replace for crosslink and protect nearest neighbors
        Args:
            model: Model object to modify
        """
        id_ = self.draw_crosslink(model)
        crosslink = model.crosslink[id_]
        
        LOG.debug(f"Processing crosslink {id_} in model {model.id}")
        LOG.debug(f"Crosslink position: {crosslink.position}")

        crosslink_state = crosslink.state
        crosslink_z = crosslink.position[2]

        if (crosslink_state not in ['replace', 'protect'] and
            self.z_min <= crosslink_z <= self.z_max):
            crosslink.state = 'replace'
            for cross in model.crosslink:
                if cross != crosslink:
                    distance = np.linalg.norm(np.array(cross.position) - np.array(crosslink.position))
                    if distance <= 10:
                        cross.state = 'replace'
                    elif 11 <= distance <= 1000:
                        cross.state = 'protect'
        
        LOG.debug(f"Finished processing crosslink {id_} in model {model.id}")

    def draw_model(self, system: Any) -> Any:
        """
        Draw a random model from system

        Args:
            system: System object to draw from

        Returns:
            Randomly selected model ID
        """
        return np.random.choice(system.get_models())

    def draw_crosslink(self, model: Any) -> int:
        """
        Draw a random crosslink from model

        Args:
            model: Model object to draw from

        Returns:
            Randomly selected crosslink index
        """
        return np.random.randint(0, len(model.crosslink))

    def write_replace(self, system: Optional[Any] = None, file: str = 'replace') -> None:
        """
        Write models with replacements to file, that is used as input for chimera

        Args:
            system: System object to write
            file: Output file name
        """
        if system is None:
            system = self.system

        with open(f"{file}.txt", 'w') as f:
            for key in system.get_models():
                for cross in system.get_model(model_id=key).crosslink:
                    if cross.state == 'replace':
                        f.write(f"{int(key)}.caps.pdb {cross.resname} {cross.resid} {cross.chain}\n")

if __name__ == "__main__":
    pass