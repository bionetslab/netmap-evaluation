from dataclasses import dataclass, field
from typing import List, Optional
import yaml
import pypiper
import sys
import os
import os.path as op
import random

sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')


from src.methods.csnet.csnet_config import CsNetConfig
from src.methods.scgenerai.scgenerai_config import ScGeneRAIConfig
from src.methods.grnboost2.grnboost2_config import GRNBoost2Config
from src.data_simulation.data_simulation_config import DataSimulationConfig
from netmap.src.utils.netmap_config import NetmapConfig
from utils import PipelineConfig
from src.utils import write_config


def read_yaml(yaml_file):
    with open(yaml_file, 'r') as f:
        data = yaml.safe_load(f)
    return data

if __name__ == "__main__":
    config_file = op.join(op.dirname(__file__), '../../configurations/pipelines/pipeline_config.yaml')
    try:
        pipeline_config = PipelineConfig.read_yaml(config_file)
    except FileNotFoundError:
        print(f"Error: {config_file} not found.")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing YAML: {e}")
        sys.exit(1)

    outfolder = op.join(pipeline_config.basedir, "benchmark_scgenerai")
    pm = pypiper.PipelineManager(name="netmap_benchmark", outfolder=outfolder)
    pm.timestamp("Hello!")



    #data_simulation_configs = read_yaml(op.join(pipeline_config.result_folder, 'all_tests.yaml'))      
    data_simulation_configs = read_yaml(op.join(pipeline_config.result_folder, 'all_tests_5_10.yaml'))      

    
    for net in data_simulation_configs:
        for netmap_trial in data_simulation_configs[net]['scgenerai']:
                for tool_trial in data_simulation_configs[net]['scgenerai'][netmap_trial]:

                    netmap_call = f"python src/methods/scgenerai/scgenerai.py --config {data_simulation_configs[net]['scgenerai'][netmap_trial][tool_trial]['config']} --dataset_config {data_simulation_configs[net]['data_simulation']['configs'][netmap_trial]}" 
                    outfile = f"{op.join(data_simulation_configs[net]['scgenerai'][netmap_trial][tool_trial]['result'], 'config.yaml')}"
                    pm.run(netmap_call, outfile )



    pm.stop_pipeline()




