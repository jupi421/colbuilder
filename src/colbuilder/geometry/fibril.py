# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import subprocess
import os
import shutil
from typing import List, Dict, Optional
from colbuilder.geometry import model

class Fibril:
    """
    Class to setup system from fibril
    """
    def __init__(self, pdb_file: Optional[str] = None, system: Optional[object] = None):
        self.pdb_file = pdb_file
        self.system = system
        self.fibril: List[str] = []
        if pdb_file:
            with open(pdb_file, 'r') as f:
                self.fibril = f.readlines()
        self.count = self.count_models()
        self.models: Dict[int, Dict] = {k: {} for k in range(self.count)}

    def seperate_system(self, pdb_file: Optional[str] = None) -> None:
        """
        Separate system in models and write pdbs
        """
        if pdb_file:
            self.pdb_file = pdb_file
            with open(pdb_file, 'r') as f:
                self.fibril = f.readlines()
        model_lines: List[str] = []
        model_cnt = 0
        last_line = ""
        for line in self.fibril:
            if line[0:3] == 'TER' and last_line[21:22] == 'C':
                self.write_model(model_cnt=model_cnt, lines=model_lines)
                model_lines = []
                model_cnt += 1
                continue
            model_lines.append(line)
            last_line = line

    def write_model(self, model_cnt: int, lines: List[str]) -> None:
        """
        Writes model to pdb-file
        """
        with open(f"{model_cnt}.caps.pdb", 'w') as f:
            for line in lines:
                f.write(line)
            f.write('END')

    def count_models(self) -> int:
        """
        Counts number of models in fibril
        """
        return sum(1 for line in self.fibril if 'TER' in line) // 3

    def build_system(self, system: Optional[object] = None) -> object:
        """
        Get models from a pdb file and checks connection
        """
        if system is None:
            system = self.system
        for model_id in range(self.count):
            model_ = model.Model(id=model_id, pdb_file=f"{model_id}.caps")
            system.add_model(model=model_)
        return system

    def write_connect(self, system: Optional[object] = None, connect_file: Optional[str] = None) -> None:
        """
        Writes system of model connections to file 
        """
        subprocess.run("rm -r N/ T/ D/ TD/ DT/", shell=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        if system is None:
            system = self.system

        with open(f"{connect_file}.txt", 'w') as f:
            for model in system.get_models():
                model_obj = system.get_model(model_id=model)
                if model_obj.connect is None:
                    continue
                elif len(model_obj.connect) == 1:
                    if not os.path.exists(os.path.join(os.getcwd(), 'N')):
                        subprocess.run("mkdir N", shell=True)
                    model_obj.type = 'N'
                    f.write(f"{model}.caps.pdb ; N\n")
                    shutil.move(os.path.join(os.getcwd(), f"{model}.caps.pdb"), os.path.join(os.getcwd(), "N"))
                else:
                    if not os.path.exists(os.path.join(os.getcwd(), model_obj.type)):
                        subprocess.run(f"mkdir {model_obj.type}", shell=True)
                    for connect in model_obj.connect:
                        f.write(f"{connect}.caps.pdb ")
                        shutil.move(os.path.join(os.getcwd(), f"{connect}.caps.pdb"), 
                                    os.path.join(os.getcwd(), model_obj.type))
                    f.write(f"; {model_obj.type}\n")

if __name__ == "__main__":
    # Add any test code or main execution here if needed
    pass