# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:43:13 2023

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


#%% Load best run for analysis - v2 prototype

data_dir = 'D:/Barbour/OneDrive - UW/Downstream/Data/20230113_Downstream/'

experiment_names = get_experiment_names(data_dir)

experiments = []

for name in experiment_names[0:7]:
    experiments.append(load_processed_data(name, data_dir))
    
    

params = {"max_diameter": 200,
          "slot_width": 150e-6 / np.sin(150. * np.pi/180) * (750/250),
          "median_particle_diameter": [10,20],
          "number of diameter bins": 51}

#%% Calculate f and process full spray


f, average_beam_length, bin_beam_length, bin_travel_time = calc_beam_length_function(experiments, params, use_scaled_diameter=False, visualize=True)
new_params = {'average beam length': average_beam_length,
              'bin travel time': bin_travel_time,
              'beam length function': f}

params.update(new_params)

pdf, VF, VFD, VFD_full, bin_volume, total_spray_area, bins, diam_dict = process_full_spray(experiments, params, use_scaled_diameter=False)

#%% Load peters data

data_dir = '/Users/mbarbour/OneDrive - UW/Downstream/Data/20230103_Downstream/'
data_dir = 'D:/Barbour/OneDrive - UW/Downstream/Data/peter_data/240SLPM_6mlmin/'

experiment_names = get_experiment_names(data_dir)

experiments_peter= []

for name in experiment_names:
    experiments_peter.append(load_processed_data(name, data_dir))

params = {"max_diameter": 200,
          "slot_width": 150e-6 / np.sin(150. * np.pi/180) * (750/250),
          "median_particle_diameter": [10,20],
          "number of diameter bins": 51}


#%% process peter's data

f, average_beam_length, bin_beam_length, bin_travel_time = calc_beam_length_function(experiments_peter, params, use_scaled_diameter=False, visualize=True)
new_params = {'average beam length': average_beam_length,
              'bin travel time': bin_travel_time,
              'beam length function': f}

params.update(new_params)

pdf_peter, VF_peter, VFD_peter, VFD_full_peter, bin_volume, total_spray_area, bins_peter, diam_dict_peter = process_full_spray(experiments_peter, params, use_scaled_diameter=False)

#%% plot the CDF of volume

fig = go.Figure()
fig.add_trace(go.Scatter(x=bins*0.467, y=np.cumsum(VFD) / np.sum(VFD),
                         marker_color=colors[0],
                         name='TV=850mL, QL=6ml/min,<br> 19 G Needle '),
              )
fig.add_trace(go.Scatter(x=bins*0.467, y=np.cumsum(VFD_peter) / np.sum(VFD_peter),
                         marker_color=colors[2],
                         name='Compressed Air, Q_air = 240 SLPM, Q_l = 6mL/min'),
              )
fig.update_yaxes(title='CDF(Volume Flux)')
fig.update_xaxes(title = 'Scaled Diameter (um)')

fig.show()


#%% How to estimate particle evaporation rate

# let's take centerline data

exp = experiments[0]
characteristics = single_measurement_statistics(exp)

d10 = characteristics['D10 scaled (um)']
U = characteristics['U mean (m/s)']

diff = 2.8e-5 # not sure sbout this either
viscosity = 15e-6 # not sure about this



Sc = viscosity / diff 
Re = U * d10 * 1e-6 / viscosity
Sh = 2 + 0.555 * Re**0.5 * Sc**(1./3)
K = 1 / (8 / 850 * np.log(1 + (.017 / (1.125 - .017)))) # need to double check

t1 = K * (d10*1e-6)**2 / diff / Sh

print("Evapoation time (s): ", t1)












