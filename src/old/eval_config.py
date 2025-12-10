from dataclasses import dataclass, field
from typing import List, Optional
import yaml

@dataclass
class PipelineConfig:
    n_top_genes: int  = 500
    clustering_method: int = 'spectral'

    @classmethod
    def read_yaml(cls, yaml_file):
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)


def write_config(c, file):
    with open(file, "w") as handle:
        yaml.safe_dump(c, handle)  
