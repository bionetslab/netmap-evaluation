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

    outfolder = op.join(pipeline_config.basedir, "benchmark_summaries_preclustered_new")
    pm = pypiper.PipelineManager(name="netmap_benchmark_summaries", outfolder=outfolder)
    pm.timestamp("Hello!")

    nr_networks = pipeline_config.number_test_modules
    data_simulation_configs = {}

    data_simulation_configs = read_yaml(op.join(pipeline_config.result_folder, 'all_tests.yaml'))      
    data_simulation_configs_2 = read_yaml(op.join(pipeline_config.result_folder, 'all_tests_5_10.yaml'))      

    data_simulation_configs = data_simulation_configs | data_simulation_configs_2

    # STEP 6. EVALUATE RESULTS
    for net in data_simulation_configs:
        # take any data simulation trial
        config_list = "--config_list "
        eval_call = f"python src/evaluation/compute_metrics.py --pipeline_config {config_file} " 

        for data_trial in data_simulation_configs[net]['data_simulation']['configs']:
            eval_call+= f"--dataset_config {data_simulation_configs[net]['data_simulation']['configs'][data_trial]} "
            for t in ['netmap', 'scgenerai', 'csnet']:
                for tool_trial in data_simulation_configs[net][t][data_trial]:
                    print(tool_trial)
                    config_list+= f"{t}_{tool_trial}={data_simulation_configs[net][t][data_trial][tool_trial]['config']} " 
            eval_call+=config_list

            #outfile = f"{op.join(pipeline_config.summary_output_dir, net, 'results.yaml')}"
            outfile = f"{op.join(pipeline_config.summary_output_dir, data_trial, net, 'clustering_score.json')}"
            print(outfile)
            pm.run(eval_call, outfile )

    for net in data_simulation_configs:
        # take any data simulation trial
        config_list = "--config_list "
        eval_call = f"python src/evaluation/compute_metrics_preclustered.py --pipeline_config {config_file} " 

        for data_trial in data_simulation_configs[net]['data_simulation']['configs']:
            eval_call+= f"--dataset_config {data_simulation_configs[net]['data_simulation']['configs'][data_trial]} "
            for t in ['grnboost2']:
                for tool_trial in data_simulation_configs[net][t][data_trial]:
                    print(tool_trial)
                    config_list+= f"{t}_{tool_trial}={data_simulation_configs[net][t][data_trial][tool_trial]['config']} " 
            eval_call+=config_list

            print(eval_call)
            outfile = f"{op.join(pipeline_config.summary_output_dir, net, 'overlaps_global.tsv')}"
            pm.run(eval_call, outfile )

    
    pm.stop_pipeline()




