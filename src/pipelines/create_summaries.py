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






def add_differential_and_common_networks(data_config, networks, basedir):
    diff_nets = networks[0:(len(networks)-1)]
    common_nets = [networks[(len(networks)-1)]]

    data_config.edgelist = [op.join( network, 'edges.tsv') for network in diff_nets]
    data_config.nodelist = [op.join( network, 'nodes.tsv') for network in diff_nets]
    
    data_config.common_edges = [op.join( network, 'edges.tsv') for network in common_nets]
    data_config.common_nodes = [op.join( network, 'nodes.tsv') for network in common_nets]

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




def pipeline_update_tool_info(config, config_path, result_path, data_trial, net, tool, tool_trial):
    if tool not in config[net]:
        config[net][tool] = {}
    if data_trial not in config[net][tool]:
        config[net][tool][data_trial] = {}

    if tool_trial not in config[net][tool][data_trial]:
        config[net][tool][data_trial][tool_trial] = {}

    config[net][tool][data_trial][tool_trial]['config']  = config_path
    config[net][tool][data_trial][tool_trial]['result']  = result_path
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

    outfolder = op.join(pipeline_config.basedir, "benchmark_summ")
    pm = pypiper.PipelineManager(name="netmap_benchmark_summaries", outfolder=outfolder)
    pm.timestamp("Hello!")

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
        networks_all = [f for f in os.listdir(pipeline_config.clustered_network_dir) if not op.isfile(op.join(pipeline_config.clustered_network_dir, f))]

        random.seed(10)
        # Generate nr_networks output data sets    
        for net in range(nr_networks):
            # For each celltype sample one network, + sample one network identical to all
            networks = random.sample(networks_all, (data_config.n_celltypes+1))

            # Wrapper function to update the data config
            data_config  = add_differential_and_common_networks(data_config, networks, basedir=pipeline_config.clustered_network_dir)
            # Name for the new data_config file
            filn = op.join(pipeline_config.simulated_data_config_dir, data_simulation_trial, f"{data_config.dataset_id}.config.yaml")
            data_config.write_yaml(yaml_file =filn)
            
            # keep track of the generated data config files in a dictionary.
            data_simulation_configs = pipeline_add_data_simulation(data_simulation_configs, 
                                                                pipeline_config.simulated_data_config_dir,
                                                                data_config.dataset_id,
                                                                data_simulation_trial)
            

    

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
                data_simulation_configs = pipeline_update_tool_info(data_simulation_configs, netmap_config_path, csnet_config.output_directory, data_trial, net, 'csnet', netmap_trial)

                csnet_config.input_data =  op.join(data_simulation_configs[net]['data_simulation']['data_path'][data_trial], 'data.h5ad')
                csnet_config.write_yaml(netmap_config_path)
            


    # STEP 5.3 NETMAP RUN
    # Create all config files for scGeneRAI run
    for net in data_simulation_configs:
        for c in pipeline_config.netmap_base_configs:
            netmap_trial = op.basename(c).replace('.yaml', '')
            netmap_config = NetmapConfig.read_yaml(yaml_file=c)

            
            for data_trial in data_simulation_configs[net]['data_simulation']['configs']:
                config_dir = op.join(pipeline_config.netmap_config_dir, netmap_trial, data_trial)
                os.makedirs(config_dir, exist_ok=True)                
                
                # add information to trial tracker
                netmap_config_path = op.join(config_dir, f"{net}.config.yaml")
                netmap_config.output_directory =  op.join(pipeline_config.result_folder, 'netmap', netmap_trial, data_trial, net)
                data_simulation_configs = pipeline_update_tool_info(data_simulation_configs, netmap_config_path, netmap_config.output_directory, data_trial, net, 'netmap', netmap_trial)

                # update config path with input data
                netmap_config.input_data =  op.join(data_simulation_configs[net]['data_simulation']['data_path'][data_trial], 'data.h5ad')
                netmap_config.write_yaml(netmap_config_path)



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
                data_simulation_configs = pipeline_update_tool_info(data_simulation_configs, netmap_config_path, scgenerai_config.output_directory, data_trial, net, 'scgenerai', netmap_trial)

                # update config path with input data
                scgenerai_config.input_data =  op.join(data_simulation_configs[net]['data_simulation']['data_path'][data_trial], 'data.h5ad')
                scgenerai_config.write_yaml(netmap_config_path)



    ## STEP 5.1 GRNBOOST2
    # Create all config files for csnet run
    for net in data_simulation_configs:
        for c in pipeline_config.grnboost2_base_configs:
            netmap_trial = op.basename(c).replace('.yaml', '')
            grnboost2_config = GRNBoost2Config.read_yaml(yaml_file=c)

            for data_trial in data_simulation_configs[net]['data_simulation']['configs']:

                config_dir = op.join(pipeline_config.grnboost2_config_dir, netmap_trial, data_trial)
                os.makedirs(config_dir, exist_ok=True)
                netmap_config_path = op.join(config_dir, f"{net}.config.yaml")
                
                grnboost2_config.output_directory =  op.join(pipeline_config.result_folder, 'grnboost2', netmap_trial, data_trial, net)
                data_simulation_configs = pipeline_update_tool_info(data_simulation_configs, netmap_config_path, grnboost2_config.output_directory, data_trial, net, 'grnboost2', netmap_trial)

                grnboost2_config.input_data =  op.join(data_simulation_configs[net]['data_simulation']['data_path'][data_trial], 'data.h5ad')
                grnboost2_config.write_yaml(netmap_config_path)



    print(data_simulation_configs)
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

            print(eval_call)
            outfile = f"{op.join(pipeline_config.summary_output_dir, net, 'results.yaml')}"
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




