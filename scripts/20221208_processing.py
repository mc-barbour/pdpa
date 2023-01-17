#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:30:33 2022

@author: mbarbour
"""

import numpy as np

import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots


from source.loader import *

pio.renderers.default = "browser"


colors = [px.colors.qualitative.Dark24[23],  px.colors.qualitative.Dark24[22],
          px.colors.qualitative.Dark2[7], px.colors.qualitative.Vivid[0]]

#%% Load experiments


data_dir = '/Users/mbarbour/OneDrive - UW/Downstream/Data/Nebulizer_v2_20221208/'

experiment_names = sorted(get_experiment_names(data_dir))

experiments = []

for name in experiment_names:

    experiments.append(load_experiment(name, data_dir))

n_exp = len(experiments)
#%% plot Diameter diff

nrows=1;ncols=n_exp
fig = make_subplots(rows=nrows, cols=ncols)

for count, exp in enumerate(experiments):

    fig.add_trace(go.Scatter(x=exp.diameter, y=exp.diameter_diff,
                             mode='markers'), row=1, col=count+1)
    fig.update_yaxes(title='Diameter difference (um)')
    fig.update_xaxes(title='Diameter (um)')

fig.show()

    
#%% scale and filter diameters - iterate and update the curve parameters in the experiment log

for exp in experiments:

    exp.filter_diameters(filter_curves=True)
    exp.scale_diameter()
    # exp.remove_zero_intensity()
    
                       
#%% plot droplet diameter vs droplet velocity

nrows=1;ncols=n_exp
fig = make_subplots(rows=nrows, cols=ncols)

for count, exp in enumerate(experiments):
    
    print(len(exp.velocity_chan1), len(exp.diameter))
    fig.add_trace(go.Scatter(x=exp.diameter, y=exp.intensity,
                             mode='markers'), row=1, col=count+1)
    fig.update_yaxes(title='Velocity (m/s)')
    fig.update_xaxes(title='Diameter (um)')

fig.show()



#%% Plot historgrams of droplet diameter
nrows=1;ncols=6
fig = make_subplots(rows=nrows, cols=ncols)

for count, exp in enumerate(experiments):
    


    fig.add_trace(go.Histogram(x=exp.diameter,
                                   xbins=dict(size=2), autobinx=False,
                                   name=exp._name,
                                   marker_color=colors[0]), row=1, col=count+1
                 )

fig.show()


#%% Compute D10

for exp in experiments:
    print("D10: {:4f}".format(np.nanmean(exp.diameter)))



#%%
def intersection(values, cutoff, BINS):

    counter = np.argmax(values >= cutoff) # find cutoff point
    point2 = np.array([BINS[counter], values[counter]]) # cutoff point
    point1 = np.array([BINS[counter-1], values[counter-1]]) # prior point
    slope = (point2[1] - point1[1])/(point2[0] - point1[0])
    intercept = point2[1] - slope*point2[0]
    dist = (cutoff - intercept) * (1/slope) # bin value exactly at cutoff - interpolated
    return dist



#%% Compute V90



nbins = 150
diameter_threshold = 150

bins = np.linspace(0,diameter_threshold,nbins)

fig = make_subplots(rows=3, cols=1)

plot_count = 0

v90 = []
measurement_distance = []
tidal_volume = []

fig = go.Figure()
for count, exp in enumerate(experiments):
    
 

    diameters = exp.diameter[~np.isnan(exp.diameter)]
    diameters = diameters[np.where(diameters < diameter_threshold)[0]]
    
    bincount = np.bincount(np.digitize(diameters,bins[0:-1]), minlength=nbins)
    
    bin_volume = 4. / 3. * np.pi * (bins / 2.) **3 * bincount # volume in each bin
    
    bin_percent_count = bincount * (1 / np.sum(bincount)) * 100
    bin_percent_volume = bin_volume * (1 / np.sum(bin_volume)) * 100
    
    bin_count_cum = np.cumsum(bin_percent_count, dtype=float)
    bin_volume_cum = np.cumsum(bin_percent_volume, dtype=float)
    
    v90.append(intersection(bin_volume_cum, 90, bins))
    measurement_distance.append(exp.measurement_distance)
    tidal_volume.append(exp.tidal_volume)
    
    fig.add_trace(go.Scatter(x=bins, y=bin_volume_cum, name=exp._name))

fig.show()
df = pd.DataFrame(data=np.array([tidal_volume, measurement_distance, v90]).T,
                  columns=['Tidal Volume', 'Measurement Distance', 'v90'])


fig = go.Figure()

fig.add_trace(go.Scatter(x=df[df['Measurement Distance']==30]['Tidal Volume'],
                         y=df[df['Measurement Distance']==30].v90,
                         name = '30mm',
                         mode = 'markers', 
                         marker_color=colors[0],
                         marker_size=20,
                         marker_symbol='circle')
                          
              )

fig.update_yaxes(title='D90 (um)')
fig.update_xaxes(title='Target tidal volume (mL)')


fig.show()


#%% 

fig = go.Figure()
fig.add_trace(go.Scatter(x=bins, y=bin_volume_cum))
fig.show()