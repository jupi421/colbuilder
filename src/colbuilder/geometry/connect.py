# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
import numpy as np
from itertools import product
from typing import Dict, List, Optional, Any
from pathlib import Path

from colbuilder.geometry import model

class Connect:
    """
    Get connections between all models in system OR between potentially added model and system.

    Attributes:
        system (Any): The system object containing models.
        pairs (Dict[float, Optional[float]]): Dictionary of model pairs.
        connect (Dict[float, List[float]]): Dictionary of connections.
        connect_file (Optional[Path]): Path to the connection file.
        external_connect (List[float]): List of external connections.
        is_line (Tuple[str, ...]): Tuple of valid line types in PDB files.
    """

    def __init__(self, system: Optional[Any] = None, connect_file: Optional[Path] = None):
        self.system = system
        self.pairs: Dict[float, Optional[float]] = {key: None for key in self.system.get_models()}
        self.connect: Dict[float, List[float]] = {}
        self.connect_file = Path(connect_file) if connect_file else None
        self.external_connect: List[float] = []
        self.is_line: Tuple[str, ...] = ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')

    def get_model_connect(self, system: Any, unit_cell: List[float]) -> bool:
        """
        Get connection between added model and all already existing models in system.

        Args:
            system (Any): The system object.
            unit_cell (List[float]): Unit cell parameters.

        Returns:
            bool: True if a connection is found, False otherwise.
        """
        transformation = system.crystal.get_t_matrix(s_matrix=unit_cell)
        add_ = model.Model(id='add', transformation=transformation, pdb_file=system.crystal.pdb_file)

        for ref_model in self.system.get_models():
            if self.get_connect(ref_model=system.get_model(model_id=ref_model), model=add_):
                return True
        return False

    def get_contact_connect(self, system: Any) -> Dict[float, List[float]]:
        """
        Get connection between all models/contacts in system.

        Args:
            system (Any): The system object.

        Returns:
            Dict[float, List[float]]: Dictionary of connections.
        """
        for ref_model, model in product(self.system.get_models(), repeat=2):
            if (ref_model != model and 
                self.get_connect(ref_model=system.get_model(model_id=ref_model),
                                 model=system.get_model(model_id=model))):
                self.pairs[ref_model] = model
        return self.merge_contacts(pairs=self.pairs)

    def merge_contacts(self, pairs: Dict[float, Optional[float]]) -> Dict[float, List[float]]:
        """
        Merges contacts to generate triplets of connections.

        Args:
            pairs (Dict[float, Optional[float]]): Dictionary of model pairs.

        Returns:
            Dict[float, List[float]]: Dictionary of merged connections.
        """
        self.connect = {key: [key] for key in pairs}

        for ref_key, key in product(pairs, repeat=2):
            if key == ref_key or pairs[key] is None or pairs[ref_key] is None:
                continue
            elif ref_key == pairs[key] or key == pairs[ref_key] or pairs[key] == pairs[ref_key]:
                self.connect[ref_key].append(key)

        return self.clean_contacts(contactpairs=self.connect)

    def get_external_connect_file(self, system: Any, connect_file: Optional[Path] = None) -> Any:
        """
        Read external connect file and update system accordingly.

        Args:
            system (Any): The system object.
            connect_file (Optional[Path]): Path to the connection file.

        Returns:
            Any: Updated system object.
        """
        if connect_file:
            self.external_connect = [float(l.split(' ')[0].replace('.caps.pdb', '')) 
                                     for l in open(connect_file.with_suffix('.txt'), 'r')]
        if np.min(self.external_connect) > 0:
            self.external_connect = [i - 1 for i in self.external_connect]

        for model_id in system.get_connect().keys():
            if model_id not in self.external_connect:
                system.get_model(model_id=model_id).connect = None
        
        return system

    def clean_contacts(self, contactpairs: Dict[float, List[float]]) -> Dict[float, List[float]]:
        """
        Clean merged contacts to prevent multi-counting in system.

        Args:
            contactpairs (Dict[float, List[float]]): Dictionary of contact pairs.

        Returns:
            Dict[float, List[float]]: Cleaned dictionary of connections.
        """
        remove_model = set([key for key in contactpairs for model in contactpairs[key] if key > model])
        for key in remove_model:
            contactpairs[key] = False
        self.connect = {k: v for k, v in contactpairs.items() if v and v is not False}
        return self.connect

    def get_connect(self, ref_model: Any, model: Any, cut_off: float = 2.0) -> bool:
        """
        Calculates distance between models: distance below cut_off (2.0 A) keep model.

        Args:
            ref_model (Any): Reference model.
            model (Any): Model to compare.
            cut_off (float): Distance cut-off in Angstroms.

        Returns:
            bool: True if models are connected, False otherwise.
        """
        for ref_c, c in product(ref_model.crosslink, model.crosslink):
            if np.linalg.norm(ref_c.position - c.position) < cut_off:
                return True
        return False

    def write_connect(self, system: Any, connect_file: Path) -> None:
        """
        Writes system of model connections to file.

        Args:
            system (Any): The system object.
            connect_file (Path): Path to the output connection file.
        """
        with open(connect_file.with_suffix('.txt'), 'w') as f:
            for model in system.get_models():
                if system.get_model(model_id=model).connect is not None:
                    if len(system.get_model(model_id=model).connect) != 1:
                        for connect in system.get_model(model_id=model).connect:
                            f.write(f"{int(connect)}.caps.pdb ")
                        f.write(f" ; {system.get_model(model_id=model).type}\n")

    def run_connect(self, system: Any, unit_cell: Optional[List[float]] = None) -> Dict[float, List[float]]:
        """
        Wrapper to determine connection between all contacts (1) or
        between an added model and the current system (2).

        Args:
            system (Any): The system object.
            unit_cell (Optional[List[float]]): Unit cell parameters.

        Returns:
            Dict[float, List[float]]: Dictionary of connections.
        """
        if unit_cell is None:
            return self.get_contact_connect(system=system)
        else:
            return self.get_model_connect(system=system, unit_cell=unit_cell)