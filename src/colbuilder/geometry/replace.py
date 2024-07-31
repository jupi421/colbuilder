# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import numpy as np
from typing import Optional, Dict, List, Any

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
            self.z_min = self.system.get_model(model_id=0.0).get_cog() - fibril_length/2
            self.z_max = self.system.get_model(model_id=0.0).get_cog() + fibril_length/2

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

        while (
            system.count_states(state='replace') / (
                system.count_states(state='none') +
                system.count_states(state='protect') +
                system.count_states(state='replace')
            ) < int(self.ratio) / 100
        ):
            model = system.get_model(model_id=self.draw_model(system=system))
            self.replace_model(model=model)
        return system

    def replace_model(self, model: Any) -> None:
        """
        Change status of replace for crosslink and protect nearest neighbors

        Args:
            model: Model object to modify
        """
        id_ = self.draw_crosslink(model)
        if (
            model.crosslink[id_].state not in ['replace', 'protect'] and
            self.z_min <= model.crosslink[id_].position[2] <= self.z_max
        ):
            model.crosslink[id_].state = 'replace'
            for cross in model.crosslink:
                if (
                    cross != model.crosslink[id_] and
                    np.linalg.norm(cross.position - model.crosslink[id_].position) <= 10
                ):
                    cross.state = 'replace'
                elif (
                    cross != model.crosslink[id_] and
                    11 <= np.linalg.norm(cross.position - model.crosslink[id_].position) <= 1000
                ):
                    cross.state = 'protect'

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