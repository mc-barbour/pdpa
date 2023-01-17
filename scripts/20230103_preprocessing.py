# -*- coding: utf-8 -*-

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

#%% Load experiments

data_dir = '/Users/mbarbour/OneDrive - UW/Downstream/Data/20230103_Downstream/'
data_dir = "../Data/20230103_Downstream/"
experiment_names = sorted(get_experiment_names(data_dir))

experiments = []

for name in experiment_names:

    experiments.append(load_raw_data(name, data_dir))

n_exp = len(experiments)

    
#%% scale and filter diameters - iterate and update the curve parameters in the experiment log

for exp in experiments[2:9]:
    exp.filter_diameters(filter_curves=True)
    #exp.scale_diameter()
    # exp.remove_zero_intensity()
    
#%% filter the particles, scale them based on surface tension, save new files

for exp in experiments:
    exp.scale_diameter()
    exp.save_processed_run()

#%% Create a few plots

fig = make_subplots(rows=1, cols=2, subplot_titles=[experiments[0]._name, experiments[10]._name])

exp_num = 0
valid_int = np.where(experiments[exp_num].valid == 1)[0]
fig.add_trace(go.Scatter(y=experiments[0].velocity_chan1[valid_int], x=experiments[0].scaled_diameter[valid_int],
                        mode='markers'), row=1, col=1)

fig.update_yaxes(title='Velocity (m/s)', row=1, col=1)
fig.update_xaxes(title='Diamter (um)', row=1, col=1)

exp_num = 10
valid_int = np.where(experiments[exp_num].valid == 1)[0]
fig.add_trace(go.Scatter(y=experiments[10].velocity_chan1[valid_int], x=experiments[10].scaled_diameter[valid_int],
                        mode='markers'), row=1, col=2)
fig.update_yaxes(title='Velocity (m/s)', row=1, col=2)
fig.update_xaxes(title='Diamter (um)', row=1, col=2)


fig.show()

#%% plot mean velocity across the jet

fig = go.Figure()
rad_pos30 = []
average_velocity30 = []

rad_pos60 = []
average_velocity60 = []

for exp in experiments:
    
    valid_int = np.where(exp.valid == 1)[0]
    if exp.centerline_distance == 30:
        rad_pos30.append(exp.radial_position)
        average_velocity30.append(np.nanmean(exp.velocity_chan1[valid_int]))
    
    if exp.centerline_distance == 60:
        rad_pos60.append(exp.radial_position)
        average_velocity60.append(np.nanmean(exp.velocity_chan1[valid_int]))

    
rad_pos60_sorted = np.array(rad_pos60)[np.argsort(rad_pos60)]
average_velocity60_sorted = np.array(average_velocity60)[np.argsort(rad_pos60)]

fig.add_trace(go.Scatter(x=rad_pos60_sorted, y=average_velocity60_sorted, name = 'y=60mm'))
fig.add_trace(go.Scatter(x=rad_pos30, y=average_velocity30, name = 'y=30mm'))
fig.update_yaxes(title='Mean Droplet Velocity (m/s)')
fig.update_xaxes(title="Radial Position (mm)")

fig.show()
        

#%% plot D10 across the jet

fig = go.Figure()
rad_pos30 = []
average_diameter30 = []

rad_pos60 = []
average_diameter60 = []

for exp in experiments:
    
    valid_int = np.where(exp.valid == 1)[0]
    if exp.centerline_distance == 30:
        rad_pos30.append(exp.radial_position)
        average_diameter30.append(np.nanmean(exp.scaled_diameter[valid_int]))
    
    if exp.centerline_distance == 60:
        rad_pos60.append(exp.radial_position)
        average_diameter60.append(np.nanmean(exp.scaled_diameter[valid_int]))

    
rad_pos60_sorted = np.array(rad_pos60)[np.argsort(rad_pos60)]
average_diameter60_sorted = np.array(average_diameter60)[np.argsort(rad_pos60)]

fig.add_trace(go.Scatter(x=rad_pos60_sorted, y=average_diameter60_sorted, name = 'y=60mm'))
fig.add_trace(go.Scatter(x=rad_pos30, y=average_diameter30, name = 'y=30mm'))
fig.update_yaxes(title='Mean Droplet Diameter, D10 (um)')
fig.update_xaxes(title="Radial Position (mm)")

fig.show()
        

 #%%
def intersection(values, cutoff, BINS):

    counter = np.argmax(values >= cutoff) # find cutoff point
    point2 = np.array([BINS[counter], values[counter]]) # cutoff point
    point1 = np.array([BINS[counter-1], values[counter-1]]) # prior point
    slope = (point2[1] - point1[1])/(point2[0] - point1[0])
    intercept = point2[1] - slope*point2[0]
    dist = (cutoff - intercept) * (1/slope) # bin value exactly at cutoff - interpolated
    return dist



#%% individual v90 plots

nbins = 150
diameter_threshold = 150

bins = np.linspace(0,diameter_threshold,nbins)

fig = make_subplots(rows=3, cols=1)

plot_count = 0

v90 = []
measurement_distance = []
tidal_volume = []

fig = go.Figure()
for count, exp in enumerate(experiments[0:10]):
    
    valid_int = np.where(exp.valid == 1)[0]
 
    diameters = exp.scaled_diameter[valid_int]
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
fig.update
fig.show()


#%% let's look at filtering sensitivity 

data_dir = '/Users/mbarbour/OneDrive - UW/Downstream/Data/20230103_Downstream/'
data_dir = "../Data/20230103_Downstream/"
experiment_names = sorted(get_experiment_names(data_dir))

experiments = []

for name in experiment_names[0:2]:

    experiments.append(load_raw_data(name, data_dir))

n_exp = len(experiments)

    
#%% scale and filter diameters - iterate and update the curve parameters in the experiment log

for exp in experiments[0:2]:
    exp.filter_diameters(filter_curves=True)
    #exp.filter_maxdiameter()
    exp.scale_diameter()
    # exp.remove_zero_intensity()

#%% compute D10 and Dv90

for exp in experiments[0:2]:

    valid_ind = np.where(exp.valid == 1)[0]
    
    diameters = exp.scaled_diameter[valid_ind]
    print("d10: {:10.6f}, dv90: {:10.6f}".format(np.mean(diameters), diameter_volume_single_meas(exp)))

#%% 


for exp in experiments[0:2]:
    valid_ind = exp.diameter < exp.maxdiameter
    
    

#     fig = go.Figure()
#     fig.add_trace(go.Scatter(x=exp.diameter, y=exp.intensity, mode='markers'))
# fig.show()






