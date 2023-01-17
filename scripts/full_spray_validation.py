# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 13:46:23 2023

@author: barbourm
"""

import numpy as np

import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots


from source.loader import *
from source.processing import *

pio.renderers.default = "browser"


colors = [px.colors.qualitative.Dark24[23],  px.colors.qualitative.Dark24[22],
          px.colors.qualitative.Dark2[7], px.colors.qualitative.Vivid[0]]


"""
Script to validate the full spray analysis code that Kee On wrote. Code conception is from Peter Huck.
Code was originally written in Matlab, and rewritten in python here
"""

#%% Load data


data_dir = '/Users/mbarbour/OneDrive - UW/Downstream/Data/keeOn_validation_set/'
data_dir = 'D:\\Barbour\\OneDrive - UW\\Downstream\\Data\\keeOn_validation_set\\'
experiment_names = sorted(get_experiment_names(data_dir))
experiment_names = get_experiment_names(data_dir)

experiments = []

for name in experiment_names:

    experiments.append(load_processed_data(name, data_dir))

n_exp = len(experiments)

exp = experiments[0]

#%% run the full spray code


params = {"max_diameter": 300,
          "slot_area": 150e-6 / np.sin(30. * np.pi/180) * (750/250),
          "median_particle_diameter": [95,115],
          "number of diameter bins": 50}


pdf, VF, VFD, VFD_full,bin_volume, total_spray_area, bins, diam_dict = process_full_spray(experiments[2::], params)


#%%
radial_pos = []

for exp in experiments:
    print(exp.radial_position)
    radial_pos.append(exp.radial_position * 1e-3)
    
dr = abs(np.mean(np.array(radial_pos[1::]) - np.array(radial_pos[0:-1])))




#%% compare pdf and vdf
pdf_matlab = np.genfromtxt(data_dir + 'PDF_code_valid.csv', delimiter=',')

fig = go.Figure()
fig.add_trace(go.Scatter(x=bins, y=pdf_matlab, name='Matlab'))
fig.add_trace(go.Scatter(x=bins,y=pdf, mode='markers', name='Python'))
fig.update_yaxes(title='PDFN(D)', type='log')
fig.update_xaxes(title='Diameter (um)')
fig.update_layout(template='plotly_white')
fig.show()


#%% vfd plot




vfd_full_matlab = np.genfromtxt(data_dir + 'VFD_code_valid.csv', delimiter=',')


VFD_matlab = np.zeros(len(bins)-1)
for bin_count in range(len(bins) - 1):
    VFD_annulus[bin_count] = np.sum(vfd_full_matlab[:, bin_count] * annulus_area)/total_spray_area



