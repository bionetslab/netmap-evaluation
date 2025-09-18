
from dataclasses import dataclass, field
from typing import List, Optional
import yaml
import os

@dataclass
class GRNBoost2Config:
    input_data: str
    tf_only: bool = False
    transcription_factors: Optional[str] = None
    grn: str = "grn.tsv"
    output_directory: str = "grnboost2_results"
    

    @classmethod
    def read_yaml(cls, yaml_file):
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)

    def write_yaml(self, yaml_file):
        with open(yaml_file, 'w') as f:
            yaml.dump(self.__dict__, f, sort_keys=False)

