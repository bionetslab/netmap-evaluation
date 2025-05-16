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
    


    data_simulation_base_configs: List[str] = field(default_factory=list)
    netmap_base_configs: List[str] = field(default_factory=list)
    csnet_base_configs: List[str] = field(default_factory=list)
    scgenerai_base_configs: List[str] = field(default_factory=list)


    @classmethod
    def read_yaml(cls, yaml_file):
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)


def write_config(c, file):
    with open(file, "w") as handle:
        yaml.safe_dump(c, handle)  




def add_differential_and_common_networks(data_config, networks):
    diff_nets = networks[0:(len(networks)-1)]
    common_nets = [networks[(len(networks)-1)]]

    data_config.edgelist = [op.join(network, 'edges.tsv') for network in diff_nets]
    data_config.nodelist = [op.join(network, 'nodes.tsv') for network in diff_nets]
    
    data_config.common_edges = [op.join(network, 'edges.tsv') for network in common_nets]
    data_config.common_nodes = [op.join(network, 'nodes.tsv') for network in common_nets]

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

    outfolder = op.join(pipeline_config.basedir, "benchmark")
    pm = pypiper.PipelineManager(name="netmap_benchmark", outfolder=outfolder)
    pm.timestamp("Hello!")


    # Step 1: pull singularity container
    target_container = "grn2gex.sif"
    singularity_pull_container = f"mkdir -p {pipeline_config.image_dir} && cd {pipeline_config.image_dir} && singularity pull {target_container} docker://hartebrodt/grn2gex && cd -"
    pm.run(cmd=singularity_pull_container, target=op.join(pipeline_config.image_dir, target_container), nofail=True)


    # Step 2: network clustering
    config_dir = op.dirname(pipeline_config.network_module_config)
    config_filename = op.basename(pipeline_config.network_module_config)
    singularity_run_command = f"mkdir -p {pipeline_config.input_network_dir} && mkdir -p {pipeline_config.clustered_network_dir} && singularity exec \
                                --bind {pipeline_config.r_script_dir}:/scripts \
                                --bind {config_dir}:/configs \
                                --bind {pipeline_config.input_network_dir}:/usr/src/app/input \
                                --bind {pipeline_config.clustered_network_dir}:/usr/src/app/output \
                                {op.join(pipeline_config.image_dir, target_container)} \
                                Rscript /scripts/network_clustering.R --config=/configs/{config_filename}"
    pm.run(singularity_run_command, pipeline_config.clustered_network_dir)

    nr_networks = pipeline_config.number_test_modules
    data_simulation_configs = {}


    #Step 3:
    #Create config files for all networks based on the base configurations
    data_simulation_configs = {}

    for c in pipeline_config.data_simulation_base_configs:

        # Read configuration, re
        data_config = DataSimulationConfig.read_yaml(yaml_file=c)
        data_simulation_trial = op.basename(c).replace('.yaml', '')
        os.makedirs(op.join(pipeline_config.simulated_data_config_dir, data_simulation_trial), exist_ok=True)

        # Get all generated clustered network files.
        networks = [f for f in os.listdir(pipeline_config.clustered_network_dir) if not op.isfile(op.join(pipeline_config.clustered_network_dir, f))]

        
        if len(networks)>nr_networks:
            random.seed(10)

            networks = random.sample(networks, nr_networks)

        # Generate nr_networks output data sets    
        for net in range(nr_networks):
            # For each celltype sample one network, + sample one network identical to all
            networks = random.sample(networks, (data_config.n_celltypes+1))

            # Wrapper function to update the data config
            data_config  = add_differential_and_common_networks(data_config, networks)
            # Name for the new data_config file
            filn = op.join(pipeline_config.simulated_data_config_dir, data_simulation_trial, f"{data_config.dataset_id}.config.yaml")
            data_config.write_yaml(yaml_file =filn)
            
            # keep track of the generated data config files in a dictionary.
            data_simulation_configs = pipeline_add_data_simulation(data_simulation_configs, 
                                                                pipeline_config.simulated_data_config_dir,
                                                                data_config.dataset_id,
                                                                data_simulation_trial)
            

        

    print(data_simulation_configs)

    # STEP 4 RUN Singularity container to generate data
    for net in data_simulation_configs:
        for trial in data_simulation_configs[net]['data_simulation']['configs']:
            trial_output_dir = op.join(pipeline_config.simulated_data_dir, trial)
            singularity_run_command = f"mkdir -p {trial_output_dir} && singularity exec \
                                            --bind {pipeline_config.r_script_dir}:/scripts \
                                            --bind {op.dirname(data_simulation_configs[net]['data_simulation']['configs'][trial])}:/configs \
                                            --bind {pipeline_config.clustered_network_dir}:/usr/src/app/input \
                                            --bind {trial_output_dir}:/usr/src/app/output \
                                            {pipeline_config.image_dir}/grn2gex.sif \
                                            Rscript /scripts/simple_data_generation.R --config=/configs/{net}.config.yaml"

            data_output_dir = op.join(trial_output_dir, str(net))
            pm.run(singularity_run_command, data_output_dir)
            
            data_simulation_configs = pipeline_update_simulation_path(data_simulation_configs, data_output_dir, trial, net)
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

                csnet_config.input_data =  op.join(data_simulation_configs[net]['data_simulation']['data_path'][data_trial], 'data.h5ad')
                csnet_config.write_yaml(netmap_config_path)
            
    for net in data_simulation_configs:
        for netmap_trial in data_simulation_configs[net]['csnet']:
                netmap_call = f"python src/methods/csnet/csnet.py --config {data_simulation_configs[net]['csnet'][netmap_trial]['config']}" 
                print(netmap_call)
                outfile = f"{op.join(data_simulation_configs[net]['csnet'][netmap_trial]['result'], 'config.yaml')}"

                pm.run(netmap_call, outfile)



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




