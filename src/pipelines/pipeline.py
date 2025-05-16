#!/usr/bin/env python

import pypiper
import sys
import os.path as op
sys.path.append(op.abspath(op.join(op.dirname(__file__), '..', '..', '..')))
import yaml
from utils.config_utils import *
import os
from pipeline_utils import *
import random

config = read_config(op.join(op.dirname(__file__), 'pipeline_config.yaml'))

outfolder = op.join(config['pipeline']['basedir'], "benchmark") # Choose a folder for your results


# Create a PipelineManager, the workhorse of pypiper
pm = pypiper.PipelineManager(name="netmap_benchmark", outfolder=outfolder)

# Timestamps to delineate pipeline sections are easy:
pm.timestamp("Hello!")


## Step 1: pull singularity container and save to the images folder
target_container = "grn2gex.sif"
singularity_pull_container = f"mydir=$(pwd) && mkdir -p {config['pipeline']['image_dir']} && cd {config['pipeline']['image_dir']} && singularity pull {target_container} docker://hartebrodt/grn2gex && cd $mydir"
pm.run(cmd = singularity_pull_container, target=target_container, nofail=True)

## Step 2: 
config_dir = op.dirname(config['pipeline']['network_module_config'])
config_file = op.basename(config['pipeline']['network_module_config'])



singularity_run_command = f"mkdir -p {config['pipeline']['input_network_dir']} && mkdir -p {config['pipeline']['clustered_network_dir']} && singularity exec \
                            --bind {config['pipeline']['r_script_dir']}:/scripts \
                            --bind {config_dir}:/configs \
                            --bind {config['pipeline']['input_network_dir']}:/usr/src/app/input \
                            --bind {config['pipeline']['clustered_network_dir']}:/usr/src/app/output \
                            {config['pipeline']['image_dir']}/grn2gex.sif \
                            Rscript /scripts/network_clustering.R --config=/configs/{config_file}"
pm.run(singularity_run_command, config['pipeline']['clustered_network_dir'])


nr_networks = config['pipeline']['number_test_modules']

#Step 3:
#Create config files for all networks based on the base configurations
data_simulation_configs = {}
for c in config['pipeline']['data_simulation_base_configs']:
    baseconfig = read_config(file=c)
    data_simulation_trial = op.basename(c).replace('.yaml', '')

    networks = [f for f in os.listdir(config['pipeline']['clustered_network_dir']) if not op.isfile(op.join(config['pipeline']['clustered_network_dir'], f))]
    os.makedirs(op.join(config['pipeline']['simulated_data_config_dir'], data_simulation_trial), exist_ok=True)

    if len(networks)>nr_networks:
        random.seed(10)
        networks = random.sample(networks, nr_networks)
    
    for net in networks:

        data_config  = add_input_to_config(baseconfig, net)
        filn = op.join(config['pipeline']['simulated_data_config_dir'], data_simulation_trial, f"{net}.config.yaml")
        os.makedirs(op.join(config['pipeline']['simulated_data_config_dir'], data_simulation_trial), exist_ok=True)

        write_config(data_config, filn)
        data_simulation_configs = pipeline_add_data_simulation(data_simulation_configs, 
                                                               config['pipeline']['simulated_data_config_dir'],
                                                               net,
                                                               data_simulation_trial)
    

print(data_simulation_configs)
## Step 4:
# Simulate data for every data simulation baseconfig



for net in data_simulation_configs:
    for trial in data_simulation_configs[net]['data_simulation']['configs']:
        trial_output_dir = op.join(config['pipeline']['simulated_data_dir'], trial)
        singularity_run_command = f"mkdir -p {trial_output_dir} && singularity exec \
                                        --bind {config['pipeline']['r_script_dir']}:/scripts \
                                        --bind {op.dirname(data_simulation_configs[net]['data_simulation']['configs'][trial])}:/configs \
                                        --bind {config['pipeline']['clustered_network_dir']}:/usr/src/app/input \
                                        --bind {trial_output_dir}:/usr/src/app/output \
                                        {config['pipeline']['image_dir']}/grn2gex.sif \
                                        Rscript /scripts/data_generation.R --config=/configs/{net}.config.yaml"

        data_output_dir = op.join(trial_output_dir, net)
        pm.run(singularity_run_command, data_output_dir)
        
        data_simulation_configs = pipeline_update_simulation_path(data_simulation_configs, data_output_dir, trial, net)
print(data_simulation_configs)




#Step 5 
# Create all config files for netmap run
for net in data_simulation_configs:
    for c in config['pipeline']['netmap_base_configs']:
        netmap_trial = op.basename(c).replace('.yaml', '')
        baseconfig = read_config(file=c)

        for data_trial in data_simulation_configs[net]['data_simulation']['configs']:
            print(data_trial)
            config_dir = op.join(config['pipeline']['netmap_config_dir'], netmap_trial, data_trial)
            os.makedirs(config_dir, exist_ok=True)
            print(config_dir)
            netmap_config_path = op.join(config_dir, f"{net}.config.yaml")
            data_simulation_configs = pipeline_update_netmap_path(data_simulation_configs, netmap_config_path, data_trial, net)
            data_simulation_configs = pipeline_update_results_path(data_simulation_configs, op.join(config['pipeline']['netmap_result_folder'], data_trial, net), data_trial, net)

            netmap_config = netmap_config_update_input_data_path(baseconfig, op.join(data_simulation_configs[net]['data_simulation']['data_path'][data_trial], 'data.h5ad'))
            netmap_config = netmap_config_update_results_data_path(baseconfig, op.join(config['pipeline']['netmap_result_folder'], netmap_trial, data_trial, net))

            write_config(netmap_config, netmap_config_path)
print(data_simulation_configs)

        
for net in data_simulation_configs:
    for netmap_trial in data_simulation_configs[net]['netmap']:
            netmap_call = f"python ../run_benchmark_20022025.py --config {data_simulation_configs[net]['netmap'][netmap_trial]} --dataset_config {data_simulation_configs[net]['data_simulation']['configs'][netmap_trial]}" 
            print(netmap_call)
            pm.run(netmap_call, f"{data_simulation_configs[net]['results'][netmap_trial]}" )


write_config(data_simulation_configs, op.join(config['pipeline']['outfolder'], 'all_tests.yaml'))

pm.stop_pipeline()