import sys
# append the path of the parent directory
sys.path.append("..")
import yaml
from utils.config_utils import *
import os.path as op
import os



def add_or_modify_constant(config, key, value):
    config['data_simulation'][key] = value
    return config


def add_input_to_config(config, current_dir):
    config['data_simulation']['edgelist'] = op.join(current_dir, 'edges.tsv')
    config['data_simulation']['nodelist'] = op.join(current_dir, 'nodes.tsv')
    config['data_simulation']['dataset_id'] = current_dir
    return config

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


def pipeline_update_netmap_path(config, config_path, trial, net):
    if 'netmap' not in config[net]:
        config[net]['netmap'] = {}

    config[net]['netmap'][trial]  = config_path
    
    return config

def pipeline_update_results_path(config, results_path, trial, net):
    if 'results' not in config[net]:
        config[net]['results'] = {}
    config[net]['results'][trial] = results_path

    return config


def netmap_config_update_input_data_path(config, input_data_path):
    config['data']['input_data'] = input_data_path
    return config
    
def netmap_config_update_results_data_path(config, results_path):
    config['results']['output_directory'] = results_path
    return config


def netmap_config_add_data_simulation(netmap_config, data_simulation_config):
    netmap_config['data_simulation'] = data_simulation_config['data_simulation']
    return netmap_config