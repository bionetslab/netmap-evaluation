from dataclasses import dataclass, field
from typing import List, Optional
import yaml
import os

@dataclass
class ScGeneRAIConfig:
    input_data: str
    tf_only: bool = False
    transcription_factors: Optional[str] = None
    adata_filename: str =  "grn_lrp.h5ad"
    grn: str = "grn_lrp.tsv"
    output_directory: str = "scgenerai_results"
    overwrite: bool  = True
    rerun: bool = True
    split: bool = False
    test_size: float = 0.3
    temp_dir: str = '/tmp'
    

    @classmethod
    def read_yaml(cls, yaml_file):
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)

    def write_yaml(self, yaml_file):
        with open(yaml_file, 'w') as f:
            yaml.dump(self.__dict__, f, sort_keys=False)

