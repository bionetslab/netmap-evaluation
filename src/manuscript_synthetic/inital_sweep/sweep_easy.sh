#!/bin/bash
pixi run python src/manuscript_synthetic/initial_sweep/parameter_sweep.py -n net_84_10865_net_88_10937_net_90_11013 -d config_easy
pixi run python src/manuscript_synthetic/initial_sweep/parameter_sweep.py -n net_133_10773_net_82_10152_net_72_10551 -d config_easy 
pixi run python src/manuscript_synthetic/initial_sweep/parameter_sweep.py -n net_53_11196_net_70_11431_net_84_9903 -d config_easy
pixi run python src/manuscript_synthetic/initial_sweep/parameter_sweep.py -n net_60_10082_net_64_11307_net_84_11226 -d config_easy
pixi run python src/manuscript_synthetic/initial_sweep/parameter_sweep.py -n net_70_11431_net_60_10082_net_135_11054 -d config_easy
pixi run python src/manuscript_synthetic/initial_sweep/parameter_sweep.py -n net_90_11013_net_75_10306_net_77_11506 -d config_easy
pixi run python src/manuscript_synthetic/initial_sweep/parameter_sweep.py -n net_172_10626_net_89_11634_net_76_10367 -d config_easy
pixi run python src/manuscript_synthetic/initial_sweep/parameter_sweep.py -n net_59_10488_net_115_11153_net_84_11226 -d config_easy  
pixi run python src/manuscript_synthetic/initial_sweep/parameter_sweep.py -n net_60_10703_net_141_11706_net_81_10586 -d config_easy  
pixi run python src/manuscript_synthetic/initial_sweep/parameter_sweep.py -n net_98_11932_net_51_10906_net_60_10082 -d config_easy
