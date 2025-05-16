from dataclasses import dataclass, field
from typing import List, Optional
import yaml


@dataclass
class CsNetConfig:
    input_data: str
    tf_only: bool = False
    transcription_factors: Optional[str] = None
    alpha: Optional[float] = None
    boxsize: Optional[float] = None
    c: Optional[float] = None
    weighted: bool = False
    filename: str = "csnet"
    format: str = ".h5ad"
    output_directory: str = "csnet_results"
    overwrite: bool = True
    normalize: bool = True

    @classmethod
    def read_yaml(cls, yaml_file):
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)

    def write_yaml(self, yaml_file):
        with open(yaml_file, 'w') as f:
            yaml.dump(self.__dict__, f, sort_keys=False)
