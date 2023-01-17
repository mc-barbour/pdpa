
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 10:39:18 2023

@author: barbourm
"""

import numpy as np
import pandas as pd

import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots


from source.loader import *
from source.processing import *

pio.renderers.default = "browser"


colors = [px.colors.qualitative.Dark24[23],  px.colors.qualitative.Dark24[22],
          px.colors.qualitative.Dark2[7], px.colors.qualitative.Vivid[0]]

#%% Load experiments

data_dir = 'D:/Barbour/OneDrive - UW/Downstream/Data/20230113_Downstream/'
experiment_names = (get_experiment_names(data_dir))

experiments = []

for name in experiment_names:

    experiments.append(load_raw_data(name, data_dir))

n_exp = len(experiments)

    
#%% scale and filter diameters - iterate and update the curve parameters in the experiment log

for exp in experiments:
    exp.filter_diameters(filter_curves=True, max_diameter=True)
    #exp.scale_diameter()
    # exp.remove_zero_intensity()
    
#%% filter the particles, scale them based on surface tension, save new files

for exp in experiments:
    exp.scale_diameter()
    exp.save_filtered_run()
    
#%% test max diameter scaling
exp = experiments[11]
valid_ind = np.where(exp.valid == 1)[0]
print(max(exp.diameter[valid_ind]))
stat_dict = single_measurement_statistics(exp, use_scaled_diameters=True)


#%% Process the experiment data and save run statistics

df_params = pd.read_excel(data_dir + 'experiment_log.xlsx', index_col='Name')
df_stats = pd.DataFrame([])
for exp in experiments:
    stat_dict = single_measurement_statistics(exp, use_scaled_diameters=True)
    stat_dict['exp name'] = exp._name
    df_run = pd.DataFrame(stat_dict, index=[0])
    df_stats = df_stats.append(df_run, ignore_index=True)

df_stats.set_index('exp name', inplace=True)

df_save = df_params.join(df_stats)

df_save.to_csv(data_dir + 'experiment_statistics.csv')
