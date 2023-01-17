#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 10:11:30 2022

@author: mbarbour
"""
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots


from source.loader import *

pio.renderers.default = "browser"


colors = [px.colors.qualitative.Dark24[23],  px.colors.qualitative.Dark24[22],
          px.colors.qualitative.Dark2[7], px.colors.qualitative.Vivid[0]]

#%%
data_dir = "/Users/mbarbour/OneDrive - UW/Downstream/Data/Neubulizer_v2_testing_20221005/"

experiment_names = get_experiment_names(data_dir)

experiment = load_experiment(experiment_names[0], data_dir)

#%%

fig = go.Figure()
fig.add_trace(go.Scatter(x=experiment.diameter, y=experiment.velocity_chan1,
              mode='markers'))
fig.show()
