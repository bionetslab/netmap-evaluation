#!/bin/bash
pixi run python src/evaluation/best_models_training.py -n net_84_10865_net_88_10937_net_90_11013 -d config_easy
pixi run python src/evaluation/best_models_training.py -n net_133_10773_net_82_10152_net_72_10551 -d config_easy 
pixi run python src/evaluation/best_models_training.py -n net_53_11196_net_70_11431_net_84_9903 -d config_easy
pixi run python src/evaluation/best_models_training.py -n net_60_10082_net_64_11307_net_84_11226 -d config_easy
pixi run python src/evaluation/best_models_training.py -n net_70_11431_net_60_10082_net_135_11054 -d config_easy
pixi run python src/evaluation/best_models_training.py -n net_90_11013_net_75_10306_net_77_11506 -d config_easy
pixi run python src/evaluation/best_models_training.py -n net_172_10626_net_89_11634_net_76_10367 -d config_easy
pixi run python src/evaluation/best_models_training.py -n net_59_10488_net_115_11153_net_84_11226 -d config_easy  
pixi run python src/evaluation/best_models_training.py -n net_60_10703_net_141_11706_net_81_10586 -d config_easy  
pixi run python src/evaluation/best_models_training.py -n net_98_11932_net_51_10906_net_60_10082 -d config_easy

pixi run python src/evaluation/best_models_training.py -n net_141_11706_net_61_10737_net_70_11670_net_51_10906 -d config_three_noise
pixi run python src/evaluation/best_models_training.py -n net_141_11706_net_81_10586_net_115_11153_net_57_10409 -d config_three_noise
pixi run python src/evaluation/best_models_training.py -n net_60_10082_net_64_11307_net_84_11226_net_133_10773 -d config_three_noise
pixi run python src/evaluation/best_models_training.py -n net_60_10082_net_90_11013_net_75_10306_net_77_11506 -d config_three_noise
pixi run python src/evaluation/best_models_training.py -n net_70_11431_net_60_10082_net_135_11054_net_53_11196 -d config_three_noise
pixi run python src/evaluation/best_models_training.py -n net_70_11431_net_84_9903_net_59_10488_net_115_11153 -d config_three_noise
pixi run python src/evaluation/best_models_training.py -n net_75_10306_net_115_11153_net_72_10551_net_113_11113 -d config_three_noise
pixi run python src/evaluation/best_models_training.py -n net_82_10152_net_72_10551_net_98_11932_net_51_10906 -d config_three_noise
pixi run python src/evaluation/best_models_training.py -n net_84_10865_net_88_10937_net_90_11013_net_60_10703 -d config_three_noise
pixi run python src/evaluation/best_models_training.py -n net_84_11226_net_172_10626_net_89_11634_net_76_10367 -d config_three_noise


pixi run python src/evaluation/best_models_training.py -n net_115_11153_net_57_10409_net_141_11706_net_61_10737_net_70_11670_net_51_10906 -d config_five
pixi run python src/evaluation/best_models_training.py -n net_70_11431_net_60_10082_net_135_11054_net_53_11196_net_74_11972_net_84_9903 -d config_five
pixi run python src/evaluation/best_models_training.py -n net_59_10488_net_115_11153_net_84_11226_net_172_10626_net_89_11634_net_76_10367 -d config_five  
pixi run python src/evaluation/best_models_training.py -n net_71_11397_net_113_11113_net_135_11054_net_53_11196_net_82_10152_net_89_11634 -d config_five
pixi run python src/evaluation/best_models_training.py -n net_60_10082_net_64_11307_net_84_11226_net_133_10773_net_82_10152_net_72_10551 -d config_five   
pixi run python src/evaluation/best_models_training.py -n net_75_10306_net_115_11153_net_72_10551_net_113_11113_net_130_11537_net_88_10937 -d config_five
pixi run python src/evaluation/best_models_training.py -n net_60_10082_net_85_11467_net_84_9903_net_72_10551_net_75_10306_net_130_10439 -d config_five
pixi run python src/evaluation/best_models_training.py -n net_84_10865_net_88_10937_net_90_11013_net_60_10703_net_141_11706_net_81_10586 -d config_five
pixi run python src/evaluation/best_models_training.py -n net_61_10737_net_153_11335_net_51_10906_net_72_10551_net_133_10773_net_70_11670 -d config_five
pixi run python src/evaluation/best_models_training.py -n net_98_11932_net_51_10906_net_60_10082_net_90_11013_net_75_10306_net_77_11506 -d config_five

