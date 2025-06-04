from dataclasses import dataclass, field
from typing import List, Optional
import yaml


@dataclass
class DataSimulationConfig:
    edgelist: List[str] = field(default_factory=list)
    nodelist: List[str] = field(default_factory=list)
    common_edges: List[str] = field(default_factory=list)
    common_nodes: List[str] = field(default_factory=list)
    n_cells: Optional[int] = None
    n_celltypes: Optional[int] = None
    seed: Optional[int] = None
    base_effect: Optional[str] = None
    mean: Optional[float] = None
    sd: Optional[float] = None
    dataset_id: Optional[str] = None
    edgelist_1: List[str] = None  
    edgelist_2: List[str] = None 
    perturbed_genes: List[str] =  None
    separator: str = '\t'
    group_key: str = 'grn'

    @classmethod
    def read_yaml(cls, yaml_file):
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)

    def write_yaml(self, yaml_file):
        with open(yaml_file, 'w') as f:
            yaml.dump(self.__dict__, f)