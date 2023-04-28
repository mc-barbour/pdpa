# -*- coding: utf-8 -*-


"""
Created on Fri Apr  7 12:43:11 2023

@author: barbourm
"""

import numpy as np

import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy.optimize import curve_fit


from source.loader import *
from source.experiment import *
from source.processing import *

pio.renderers.default = "browser"


colors = [px.colors.qualitative.Dark24[23],  px.colors.qualitative.Dark24[22],
          px.colors.qualitative.Dark2[7], px.colors.qualitative.Vivid[0]]


#%% Load experiments

data_dir = 'D:/Barbour/OneDrive - UW/Downstream/Data/20230412_Downstream/'
experiment_names = (get_experiment_names(data_dir))

experiments = []

for name in experiment_names:

    experiments.append(load_raw_data(name, data_dir))

n_exp = len(experiments)


#%% Set the filter curves
for exp in experiments:
    exp.filter_diameters(filter_curves=True, max_diameter=True, show_plot=True)
    exp.estimate_flowrate()
    exp.save_filtered_run()
    
    
#%% Calculate per run statistics and save
df_params = pd.read_excel(data_dir + 'experiment_log.xlsx', index_col='Name')
df_stats = pd.DataFrame([])
for exp in experiments:
    stat_dict = single_measurement_statistics(exp, use_scaled_diameters=False)
    stat_dict['exp name'] = exp._name
    df_run = pd.DataFrame(stat_dict, index=[0])
    df_stats = df_stats.append(df_run, ignore_index=True)

df_stats.set_index('exp name', inplace=True)

df_save = df_params.join(df_stats)

df_save.to_csv(data_dir + 'experiment_statistics.csv')


#%% Generate plots - load the downstream data

df_tmp = pd.read_csv('D:/Barbour/OneDrive - UW/Downstream/Data/20230405_Downstream/experiment_statistics.csv')

df_ethanol_0405 = df_tmp[df_tmp.Liquid=='ethanol']


df_ethanol_0412 = df_save.iloc[0:4]
fig = make_subplots(rows=1, cols=3)

fig.add_trace(go.Scatter(x=df_ethanol_0412['Radial position (mm)'], y=df_ethanol_0412['U mean (m/s)'],
                         marker_color=colors[2], name='40 cm tubing'),row=1,col=1)
fig.add_trace(go.Scatter(x=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==30]['Radial position (mm)'],
                         y=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==30]['U mean (m/s)'],
                         marker_color=colors[0], name='30 mm open air'),row=1,col=1)
fig.add_trace(go.Scatter(x=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==90]['Radial position (mm)'], 
                         y=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==90]['U mean (m/s)'],
                         marker_color=colors[1], name='90 mm open air'),row=1,col=1)
fig.update_yaxes(title="Velocity (m/s)", row=1, col=1)





fig.add_trace(go.Scatter(x=df_ethanol_0412['Radial position (mm)'], y=df_ethanol_0412['D10 (um)'],
                         marker_color=colors[2],showlegend=False),row=1,col=2)

fig.add_trace(go.Scatter(x=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==30]['Radial position (mm)'],
                         y=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==30]['D10 (um)'],
                         marker_color=colors[0],showlegend=False),row=1,col=2)
fig.add_trace(go.Scatter(x=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==90]['Radial position (mm)'], 
                         y=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==90]['D10 (um)'],
                         marker_color=colors[1],showlegend=False),row=1,col=2)
fig.update_yaxes(title="D10 (um)", row=1, col=2)




fig.add_trace(go.Scatter(x=df_ethanol_0412['Radial position (mm)'], y=df_ethanol_0412['D32 (um)'],
                         marker_color=colors[2],showlegend=False),row=1,col=3)

fig.add_trace(go.Scatter(x=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==30]['Radial position (mm)'],
                         y=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==30]['D32 (um)'],
                         marker_color=colors[0],showlegend=False),row=1,col=3)
fig.add_trace(go.Scatter(x=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==90]['Radial position (mm)'], 
                         y=df_ethanol_0405[df_ethanol_0405["Centerline distance (mm)"]==90]['D32 (um)'],
                         marker_color=colors[1],showlegend=False),row=1,col=3)
fig.update_yaxes(title="D32 (um)", row=1, col=3)

fig.update_xaxes(title='Radial Position (mm)')



fig.update_xaxes(title='Radial Position (mm)')


fig.show()


#%% now, let's look at burst data - going to load the filtered data

experiment_names = get_experiment_names(data_dir)

experiments = []

for name in experiment_names:
    experiments.append(load_processed_data(name, data_dir))
#%%

exp = experiments[0]

valid_ind = np.where(exp.valid == 1)[0]

fig = go.Figure()
fig.add_trace(go.Scatter(y=exp.velocity_chan1[valid_ind],
                         x=exp.time_chan1[valid_ind],
                         mode='markers'))

fig.show()



#%% define burst windows
min_sep = 1.0

time = exp.time_chan1

burst_ind = np.where(np.diff(time) > min_sep)[0]

burst_length = []
for count in range(len(burst_ind)-1):
    dt = time[burst_ind[count+1]] - time[burst_ind[count]+1]

    burst_length.append(dt)


#%% Lets interpolate onto a single burst sequence and then average.



