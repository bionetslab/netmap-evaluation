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
from src.data_simulation.data_simulation_config import DataSimulationConfig
from netmap.src.utils.netmap_config import NetmapConfig

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
    scgenerai_config_dir: str
    perturb_seq_config_dir: str
    perturb_seq_subset_dir: str



    data_simulation_base_configs: List[str] = field(default_factory=list)
    netmap_base_configs: List[str] = field(default_factory=list)
    csnet_base_configs: List[str] = field(default_factory=list)
    scgenerai_base_configs: List[str] = field(default_factory=list)
    perturb_seq_base_configs: List[str] = field(default_factory=list)


    @classmethod
    def read_yaml(cls, yaml_file):
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)


def write_config(c, file):
    with open(file, "w") as handle:
        yaml.safe_dump(c, handle)  




def add_differential_and_common_networks(data_config, networks, basedir):
    diff_nets = networks[0:(len(networks)-1)]
    common_nets = [networks[(len(networks)-1)]]

    data_config.edgelist = [op.join(basedir, network, 'edges.tsv') for network in diff_nets]
    data_config.nodelist = [op.join(basedir, network, 'nodes.tsv') for network in diff_nets]
    
    data_config.common_edges = [op.join(basedir, network, 'edges.tsv') for network in common_nets]
    data_config.common_nodes = [op.join(basedir, network, 'nodes.tsv') for network in common_nets]

    data_config.dataset_id  = '_'.join(networks)

    
    return data_config


def pipeline_add_data_simulation(configuration_files, directory, net, trial):
    filn  = op.join(directory, trial, f"{net}.config.yaml")
    if not net in configuration_files:
        configuration_files[net] = {}
    if not 'data_simulation' in configuration_files[net]:
        configuration_files[net]['data_simulation'] = {}
    if not 'configs' in configuration_files[net]['data_simulation']:
        configuration_files[net]['data_simulation']['configs'] = {trial: filn}
    else:
        configuration_files[net]['data_simulation']['configs'][trial] = filn
    return configuration_files

def pipeline_add_perturb_seq(configuration_files, directory, net, trial):
    filn  = op.join(directory, f"{net}.config.yaml")
    if not net in configuration_files:
        configuration_files[net] = {}
    if not 'data_simulation' in configuration_files[net]:
        configuration_files[net]['data_simulation'] = {}
    if not 'configs' in configuration_files[net]['data_simulation']:
        configuration_files[net]['data_simulation']['configs'] = {trial: filn}
    else:
        configuration_files[net]['data_simulation']['configs'][trial] = filn
    return configuration_files

def pipeline_update_simulation_path(config, trial_output_dir, trial, net):
    if 'data_path' not in config[net]['data_simulation']:
        config[net]['data_simulation']['data_path'] = {trial: trial_output_dir}
    else:
        config[net]['data_simulation']['data_path'][trial] = trial_output_dir
    return config




def pipeline_update_tool_info(config, config_path, result_path, trial, net, tool):
    if tool not in config[net]:
        config[net][tool] = {}
    if trial not in config[net][tool]:
        config[net][tool][trial] = {}


    config[net][tool][trial]['config']  = config_path
    config[net][tool][trial]['result']  = result_path
    return config




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

    outfolder = op.join(pipeline_config.basedir, "benchmark_perturbation")
    pm = pypiper.PipelineManager(name="netmap_benchmark", outfolder=outfolder)
    pm.timestamp("Hello!")


    data_simulation_configs = {}
    #STEP 4a: Generate Data from perturb seq study
    data_generation_run_command = "python src/data_simulation/create_perturbseq_data.py --config configurations/perturb_seq/config.yaml"
    pm.run(data_generation_run_command, op.join(pipeline_config.result_folder, "perturb_seq/clustering_metrics_2tfs.tsv"))

    perturb_seq_configurations = os.listdir(op.join(pipeline_config.result_folder, 'configurations',  "perturb_seq"))
    print(perturb_seq_configurations)
    
    for ps in perturb_seq_configurations:
        config_basename =op.basename(ps).replace('.config.yaml', '')
        data_simulation_configs = pipeline_add_perturb_seq(data_simulation_configs, pipeline_config.perturb_seq_config_dir, config_basename, "perturb_seq")
        print(data_simulation_configs)
        data_simulation_configs = pipeline_update_simulation_path(data_simulation_configs, op.join(pipeline_config.perturb_seq_subset_dir, config_basename), 'perturb_seq', config_basename)

    print(data_simulation_configs)



    # STEP 5: generate tool configs.

    ## STEP 5.1 CSNET
    # Create all config files for csnet run
    for net in data_simulation_configs:
        for c in pipeline_config.csnet_base_configs:
            netmap_trial = op.basename(c).replace('.yaml', '')
            csnet_config = CsNetConfig.read_yaml(yaml_file=c)

            for data_trial in data_simulation_configs[net]['data_simulation']['configs']:

                config_dir = op.join(pipeline_config.csnet_config_dir, netmap_trial, data_trial)
                os.makedirs(config_dir, exist_ok=True)
                netmap_config_path = op.join(config_dir, f"{net}.config.yaml")
                
                csnet_config.output_directory =  op.join(pipeline_config.result_folder, 'csnet', netmap_trial, data_trial, net)
                data_simulation_configs = pipeline_update_tool_info(data_simulation_configs, netmap_config_path, csnet_config.output_directory, data_trial, net, 'csnet')

                print(op.join(data_simulation_configs[net]['data_simulation']['data_path'][data_trial], 'data.h5ad'))
                csnet_config.input_data =  op.join(data_simulation_configs[net]['data_simulation']['data_path'][data_trial], 'data.h5ad')
                csnet_config.write_yaml(netmap_config_path)
            
    for net in data_simulation_configs:
        for netmap_trial in data_simulation_configs[net]['csnet']:
                netmap_call = f"python src/methods/csnet/csnet.py --config {data_simulation_configs[net]['csnet'][netmap_trial]['config']}" 
                print(netmap_call)
                outfile = f"{op.join(data_simulation_configs[net]['csnet'][netmap_trial]['result'], 'config.yaml')}"

                pm.run(netmap_call, outfile)



    # STEP 5.3 NETMAP RUN
    # Create all config files for scGeneRAI run
    for net in data_simulation_configs:
        for c in pipeline_config.netmap_base_configs:
            netmap_trial = op.basename(c).replace('.yaml', '')
            print(c)
            netmap_config = NetmapConfig.read_yaml(yaml_file=c)

            for data_trial in data_simulation_configs[net]['data_simulation']['configs']:
                config_dir = op.join(pipeline_config.netmap_config_dir, netmap_trial, data_trial)
                os.makedirs(config_dir, exist_ok=True)                
                
                # add information to trial tracker
                netmap_config_path = op.join(config_dir, f"{net}.config.yaml")
                netmap_config.output_directory =  op.join(pipeline_config.result_folder, 'netmap', netmap_trial, data_trial, net)
                data_simulation_configs = pipeline_update_tool_info(data_simulation_configs, netmap_config_path, netmap_config.output_directory, data_trial, net, 'netmap')

                # update config path with input data
                netmap_config.input_data =  op.join(data_simulation_configs[net]['data_simulation']['data_path'][data_trial], 'data.h5ad')
                netmap_config.write_yaml(netmap_config_path)


    for net in data_simulation_configs:
        for netmap_trial in data_simulation_configs[net]['netmap']:
                netmap_call = f"python src/methods/netmap/netmap_runner.py --config {data_simulation_configs[net]['netmap'][netmap_trial]['config']} --dataset_config {data_simulation_configs[net]['data_simulation']['configs'][netmap_trial]}" 
                print(netmap_call)
                outfile = f"{op.join(data_simulation_configs[net]['netmap'][netmap_trial]['result'], 'config.yaml')}"
                pm.run(netmap_call, outfile )


    # STEP 5.3 SCGENERAI RUN
    # Create all config files for scGeneRAI run
    for net in data_simulation_configs:
        for c in pipeline_config.scgenerai_base_configs:
            netmap_trial = op.basename(c).replace('.yaml', '')
            scgenerai_config = ScGeneRAIConfig.read_yaml(yaml_file=c)

            for data_trial in data_simulation_configs[net]['data_simulation']['configs']:
                config_dir = op.join(pipeline_config.scgenerai_config_dir, netmap_trial, data_trial)
                os.makedirs(config_dir, exist_ok=True)                
                
                # add information to trial tracker
                netmap_config_path = op.join(config_dir, f"{net}.config.yaml")
                scgenerai_config.output_directory =  op.join(pipeline_config.result_folder, 'scgenerai', netmap_trial, data_trial, net)
                data_simulation_configs = pipeline_update_tool_info(data_simulation_configs, netmap_config_path, scgenerai_config.output_directory, data_trial, net, 'scgenerai')

                # update config path with input data
                scgenerai_config.input_data =  op.join(data_simulation_configs[net]['data_simulation']['data_path'][data_trial], 'data.h5ad')
                scgenerai_config.write_yaml(netmap_config_path)


    for net in data_simulation_configs:
        for netmap_trial in data_simulation_configs[net]['scgenerai']:
                netmap_call = f"python src/methods/scgenerai/scgenerai.py --config {data_simulation_configs[net]['scgenerai'][netmap_trial]['config']} --dataset_config {data_simulation_configs[net]['data_simulation']['configs'][netmap_trial]}" 
                print(netmap_call)
                outfile = f"{op.join(data_simulation_configs[net]['scgenerai'][netmap_trial]['result'], 'config.yaml')}"
                pm.run(netmap_call, outfile )



    write_config(data_simulation_configs, file=op.join(pipeline_config.outfolder, 'all_tests.yaml'))

    pm.stop_pipeline()




