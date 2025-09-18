import os
import sys
os.environ["CUDA_VISIBLE_DEVICES"]="0"


import torch
torch.cuda.is_available()

## External imports
import os
import pandas as pd
import numpy as np
import scanpy as sc
import time
import torch

import warnings
warnings.filterwarnings("ignore")
import os.path as op



import pandas as pd

sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')
from netmap.src.utils.data_utils import *
from netmap.src.utils.tf_utils import *
from netmap.src.utils.netmap_config import NetmapConfig
from netmap.src.model.nbautoencoder import *
from netmap.src.model.nbautoencoder import train_autoencoder
from netmap.src.old.inferrence_simple import *
from netmap.src.old.model_concept_v2 import *

from src.data_simulation.data_simulation_config import DataSimulationConfig




def parse_args():
    parser = argparse.ArgumentParser(description="Argument parser with two configuration files.")
    
    parser.add_argument(
        '--config',
        type=str,
        default='configuration.json',
        help='Path to the main configuration file (default: configuration.json)'
    )
    parser.add_argument(
        '--dataset_config',
        type=str,
        default='dataset_configuration.json',
        help='Path to the dataset configuration file (default: dataset_configuration.json)'
    )
    
    return parser.parse_args()

if __name__ == "__main__":
    import argparse

    args = parse_args()
    print(f"Main Configuration File: {args.config}")
    print(f"Dataset Configuration File: {args.dataset_config}")

    
    config = NetmapConfig.read_yaml(args.config)
    dataset_config = DataSimulationConfig.read_yaml(args.dataset_config)

    
    run_netmap(config, dataset_config)