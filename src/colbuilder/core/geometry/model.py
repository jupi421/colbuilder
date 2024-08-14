# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations
import numpy as np
from typing import List, Optional, Any

from colbuilder.core.geometry import crosslink

class Model:
    """
    Class for model that stores all relevant information about it.

    Attributes:
        id (Any): Identifier for the model.
        transformation (List[float]): Transformation matrix.
        unit_cell (Optional[List[float]]): Unit cell parameters.
        connect (Optional[List[float]]): List of connected model IDs.
        connect_id (Optional[float]): ID of the connection.
        crosslink (List[crosslink.Crosslink]): List of crosslink objects.
        type (str): Type of the model based on crosslinks.
        cog (np.ndarray): Center of geometry of the model.
    """

    def __init__(self, id: Any, transformation: List[float], unit_cell: Optional[List[float]] = None,
                 connect: Optional[List[float]] = None, connect_id: Optional[float] = None,
                 pdb_file: Optional[str] = None):
        self.id = id
        self.transformation = transformation
        self.unit_cell = unit_cell
        self.connect = connect
        self.connect_id = connect_id
        self.crosslink = self.add_crosslink(crosslink=crosslink.read_crosslink(pdb_file=pdb_file))
        self.type = "".join(set(cross.type for cross in self.crosslink))
        self.cog = self.get_cog()

    def add_connect(self, connect_id: float, connect: List[float]) -> None:
        """
        Add information about model's connections.

        Args:
            connect_id (float): ID of the connection.
            connect (List[float]): List of connected model IDs.
        """
        self.connect_id = connect_id
        self.connect = connect

    def delete_connect(self, connect_id: float) -> None:
        """
        Delete connect_id from model in system.

        Args:
            connect_id (float): ID of the connection to delete.
        """
        if self.connect:
            self.connect = [i for i in self.connect if i != connect_id]

    def add_crosslink(self, crosslink: List[crosslink.Crosslink]) -> List[crosslink.Crosslink]:
        """
        Add and transform crosslink according to transformation matrix of model.

        Args:
            crosslink (List[crosslink.Crosslink]): List of crosslink objects.

        Returns:
            List[crosslink.Crosslink]: Transformed list of crosslink objects.
        """
        for cross in crosslink:
            cross.set_transform(model_id=self.id, transform=self.transformation)
        return crosslink

    def get_cog(self) -> np.ndarray:
        """
        Get center-of-geometry of model.

        Returns:
            np.ndarray: Center of geometry.
        """
        return np.mean([cross.position for cross in self.crosslink], axis=0)

    def count_state(self, state: str) -> int:
        """
        Count crosslinks with certain state (no, replace, protect) in model.

        Args:
            state (str): State to count ('no', 'replace', or 'protect').

        Returns:
            int: Number of crosslinks with the specified state.
        """
        return sum(1 for cross in self.crosslink if cross.state == state)