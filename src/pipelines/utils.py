from dataclasses import dataclass, field
from typing import List, Optional
import yaml

@dataclass
class PipelineConfig:
    basedir: str
    outfolder: str
    image_dir: str
    r_script_dir: str
    network_module_config: str
    number_test_modules: int
    input_network_dir: str
    clustered_network_dir: str
    simulated_data_dir: str
    simulated_data_config_dir: str
    result_folder: str
    netmap_config_dir: str
    csnet_config_dir: str
    grnboost2_config_dir: str
    scgenerai_config_dir: str
    perturb_seq_config_dir: str
    perturb_seq_subset_dir: str
    summary_output_dir: str



    data_simulation_base_configs: List[str] = field(default_factory=list)
    netmap_base_configs: List[str] = field(default_factory=list)
    csnet_base_configs: List[str] = field(default_factory=list)
    scgenerai_base_configs: List[str] = field(default_factory=list)
    grnboost2_base_configs: List[str] = field(default_factory=list)

    perturb_seq_base_configs: List[str] = field(default_factory=list)



    @classmethod
    def read_yaml(cls, yaml_file):
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)


def write_config(c, file):
    with open(file, "w") as handle:
        yaml.safe_dump(c, handle)  
