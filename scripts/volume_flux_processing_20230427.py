# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 10:55:19 2023

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


#%% First let's process the volume flux data from 0412

data_dir = 'D:/Barbour/OneDrive - UW/Downstream/Data/20230412_Downstream/'

experiment_names = get_experiment_names(data_dir)

experiments = []

for name in experiment_names:
    experiments.append(load_processed_data(name, data_dir))
    

#%%

params = {"max diameter": 54,
          "slot width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "number of diameter bins": 51,
          "beam length cutoff": 95,
          "function type": 'plateau', 
          "spray mode": 'burst',
          "minimum burst seperation": 1.0}

beam_length_func_params = calc_diameter_beam_length(experiments[0:4], params, visualize=True)

pdfN, volume_flux_spray, volume_flux_position, bin_volume, total_spray_area, bins, diam_dict = full_spray_metrics(experiments[0:4], beam_length_func_params, params)

ethanol_tubing_dict = {"pdf": pdfN,
                       'volume flux': volume_flux_spray,
                       'volume flux position':  volume_flux_position,
                       'bin volume':bin_volume,
                       'total measured area': total_spray_area,
                       'bins': bins}

ethanol_tubing_dict.update(diam_dict)

fig = go.Figure()
fig.add_trace(go.Scatter(y=np.sum(volume_flux_position, axis=1)))
fig.update_yaxes(title='Volume flux (m3 / (m2 s))')
fig.show()


#%% Let's look at the ethanol data from 0405
data_dir = 'D:/Barbour/OneDrive - UW/Downstream/Data/20230405_Downstream/'

experiment_names = get_experiment_names(data_dir)

experiments = []

for count, name in enumerate(experiment_names):
    print(count, name)
    experiments.append(load_processed_data(name, data_dir))

#extract the ethanol data
exp_ethanol_30mm = experiments[13:18]
exp_ethanol_90mm = experiments[18:25]

exp_hexane_30mm = experiments[25:30]
exp_hexane_90mm = experiments[30:36]



#%%

for exp in exp_hexane_30mm:
    print(exp.radial_position)
for exp in exp_hexane_90mm:
    print(exp.radial_position)

#%% process the 30mm ethanol data


params = {"max diameter": 90,
          "slot width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "number of diameter bins": 51,
          "beam length cutoff": 98,
          "function type": 'plateau', 
          "spray mode": 'burst',
          "minimum burst seperation": 1.0}

beam_length_func_params = calc_diameter_beam_length(exp_ethanol_30mm, params, visualize=True)

pdfN, volume_flux_spray, volume_flux_position, bin_volume, total_spray_area, bins, diam_dict = full_spray_metrics(exp_ethanol_30mm, beam_length_func_params, params)

ethanol_30mm_dict = {"pdf": pdfN,
                       'volume flux': volume_flux_spray,
                       'volume flux position':  volume_flux_position,
                       'bin volume':bin_volume,
                       'total measured area': total_spray_area,
                       'bins': bins}

ethanol_30mm_dict.update(diam_dict)

fig = go.Figure()
fig.add_trace(go.Scatter(y=np.sum(volume_flux_position, axis=1)))
fig.update_yaxes(title='Volume flux (m3 / (m2 s))')
fig.show()


#%% the center measurement has only one measurement - because of the pdpa settings

params = {"max diameter": 90,
          "slot width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "number of diameter bins": 51,
          "beam length cutoff": 95,
          "function type": 'plateau', 
          "spray mode": 'continuous',
          "minimum burst seperation": 1.0}

pdfN, volume_flux_spray, volume_flux_position0, bin_volume, total_spray_area, bins, diam_dict = full_spray_metrics(exp_ethanol_90mm[0:2], beam_length_func_params, params)
vol_flux0 = np.sum(volume_flux_position0, axis=1)[0]


#%% process the 90mm ethanol data

params = {"max diameter": 90,
          "slot width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "number of diameter bins": 51,
          "beam length cutoff": 95,
          "function type": 'plateau', 
          "spray mode": 'burst',
          "minimum burst seperation": 1.0}

beam_length_func_params = calc_diameter_beam_length(exp_ethanol_90mm, params, visualize=True)

pdfN, volume_flux_spray, volume_flux_position, bin_volume, total_spray_area, bins, diam_dict = full_spray_metrics(exp_ethanol_90mm[1::], beam_length_func_params, params)

ethanol_90mm_dict = {"pdf": pdfN,
                       'volume flux': volume_flux_spray,
                       'volume flux position':  volume_flux_position,
                       'bin volume':bin_volume,
                       'total measured area': total_spray_area,
                       'bins': bins}

ethanol_90mm_dict.update(diam_dict)

fig = go.Figure()
fig.add_trace(go.Scatter(y=np.sum(volume_flux_position, axis=1)))
fig.update_yaxes(title='Volume flux (m3 / (m2 s))')
fig.show()

#%% 
exp = exp_ethanol_90mm[0]
exp._name
valid_int = np.where(exp.valid == 1)[0]
time = exp.time_chan1[valid_int]
d = exp.diameter[valid_int]

fig = go.Figure()
fig.add_trace(go.Scatter(x=time, y=d))
fig.show()



#%% Plot the radial profiles of volume flux

radial_pos_tubing = np.array([0,2,4,6])
radial_pos_30mm = np.array([0,2,4,6,8])
radial_pos_90mm = np.array([0,2,4,6,8,10,12])


vol_flux_radial_90 = np.sum(ethanol_90mm_dict['volume flux position'], axis=1)
vol_flux_radial_90 = np.insert(vol_flux_radial_90, 0, vol_flux0)


fig = go.Figure()
fig.add_trace(go.Scatter(x=radial_pos_tubing, y=np.sum(ethanol_tubing_dict['volume flux position'], axis=1),
                         marker_color=colors[0],
                         name='2 Sections Tubing'))
fig.add_trace(go.Scatter(x=radial_pos_30mm, y=np.sum(ethanol_30mm_dict['volume flux position'], axis=1),
              marker_color=colors[1],
              name='30mm Open Air'))
fig.add_trace(go.Scatter(x=radial_pos_90mm, y=vol_flux_radial_90,
              marker_color=colors[2],
              name='90mm Open Air'))

fig.update_yaxes(title='Volume flux (m3 / (m2 s))')
fig.update_xaxes(title='Radial Position (mm)')

fig.show()

print("Spray Volume @ 30mm: {:.4f} ml/min".format(np.sum(ethanol_30mm_dict['volume flux']) *60e6))
print("Spray Volume @ 90mm: {:.4f} ml/min".format(np.sum(ethanol_90mm_dict['volume flux']) *60e6))

print("Spray Volume @ tubing: {:.4f} ml/min".format(np.sum(ethanol_tubing_dict['volume flux']) *60e6))


#%% lets' try and estimate the liquid volume flux of ethanol

spray_volume_rate_30mm = np.sum(ethanol_30mm_dict['volume flux']) *60e6



radial_pos, annulus_area, total_spray_area = compute_annulus_area(exp_ethanol_90mm)
spray_volume_rate_90mm = np.sum(vol_flux_radial_90*annulus_area)*60e6

spray_volume_rate_tubing = np.sum(ethanol_tubing_dict['volume flux']) *60e6



#%% now lets look at hexane - 30mm

params = {"max diameter": 78,
          "slot width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "number of diameter bins": 51,
          "beam length cutoff": 95,
          "function type": 'plateau', 
          "spray mode": 'burst',
          "minimum burst seperation": 1.0}

beam_length_func_params = calc_diameter_beam_length(exp_hexane_30mm, params, visualize=True)

pdfN, volume_flux_spray, volume_flux_position, bin_volume, total_spray_area, bins, diam_dict = full_spray_metrics(exp_hexane_30mm, beam_length_func_params, params)

hexane_30mm_dict = {"pdf": pdfN,
                       'volume flux': volume_flux_spray,
                       'volume flux position':  volume_flux_position,
                       'bin volume':bin_volume,
                       'total measured area': total_spray_area,
                       'bins': bins}

hexane_30mm_dict.update(diam_dict)

fig = go.Figure()
fig.add_trace(go.Scatter(y=np.sum(volume_flux_position, axis=1)))
fig.update_yaxes(title='Volume flux (m3 / (m2 s))')
fig.show()



#%% now lets look at hexane - 90mm

params = {"max diameter": 78,
          "slot width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "number of diameter bins": 51,
          "beam length cutoff": 95,
          "function type": 'plateau', 
          "spray mode": 'burst',
          "minimum burst seperation": 1.0}

beam_length_func_params = calc_diameter_beam_length(exp_hexane_90mm, params, visualize=True)

pdfN, volume_flux_spray, volume_flux_position, bin_volume, total_spray_area, bins, diam_dict = full_spray_metrics(exp_hexane_90mm, beam_length_func_params, params)

hexane_90mm_dict = {"pdf": pdfN,
                       'volume flux': volume_flux_spray,
                       'volume flux position':  volume_flux_position,
                       'bin volume':bin_volume,
                       'total measured area': total_spray_area,
                       'bins': bins}

hexane_90mm_dict.update(diam_dict)

fig = go.Figure()
fig.add_trace(go.Scatter(y=np.sum(volume_flux_position, axis=1)))
fig.update_yaxes(title='Volume flux (m3 / (m2 s))')
fig.show()


#%% plot hexane decay

radial_pos_30mm = np.array([0,2,4,6,8])
radial_pos_90mm = np.array([0,2,4,6,8,10])


fig = go.Figure()

fig.add_trace(go.Scatter(x=radial_pos_30mm, y=np.sum(hexane_30mm_dict['volume flux position'], axis=1),
              marker_color=colors[1],
              name='30mm Open Air'))
fig.add_trace(go.Scatter(x=radial_pos_90mm, y=np.sum(hexane_90mm_dict['volume flux position'], axis=1),
              marker_color=colors[2],
              name='90mm Open Air'))

fig.update_yaxes(title='Volume flux (m3 / (m2 s))')
fig.update_xaxes(title='Radial Position (mm)')

fig.show()

print("Spray Volume @ 30mm: {:.4f}".format(np.sum(hexane_30mm_dict['volume flux']) *60e6))
print("Spray Volume @ 90mm: {:.4f}".format(np.sum(hexane_90mm_dict['volume flux']) *60e6))



#%% Let's load the data from 4/24

data_dir = 'D:/Barbour/OneDrive - UW/Downstream/Data/20230424_Downstream/'

experiment_names = get_experiment_names(data_dir)

experiments = []

for count, name in enumerate(experiment_names):
    print(count, name)
    experiments.append(load_processed_data(name, data_dir))

#%%
params = {"max diameter": 50,
          "slot width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "number of diameter bins": 51,
          "beam length cutoff": 95,
          "function type": 'plateau', 
          "spray mode": 'burst',
          "minimum burst seperation": 3.0}

beam_length_func_params = calc_diameter_beam_length(experiments, params, visualize=True)

pdfN, volume_flux_spray, volume_flux_position, bin_volume, total_spray_area, bins, diam_dict = full_spray_metrics(experiments[0:4], beam_length_func_params, params)




random_dict = {"pdf": pdfN,
                       'volume flux': volume_flux_spray,
                       'volume flux position':  volume_flux_position,
                       'bin volume':bin_volume,
                       'total measured area': total_spray_area,
                       'bins': bins}
random_dict.update(diam_dict)



fig = go.Figure()
fig.add_trace(go.Scatter(y=np.sum(volume_flux_position, axis=1)))
fig.update_yaxes(title='Volume flux (m3 / (m2 s))')
fig.show()




#%% debug calc-time


exp = experiments[3]
exp._name
valid_int = np.where(exp.valid == 1)[0]
time = exp.time_chan1[valid_int]
d = exp.diameter[valid_int]

fig = go.Figure()
fig.add_trace(go.Scatter(x=time, y=d))
fig.show()


#%% estimate the liquid volume flux using a linear decay

n_pos = 10
radial_pos = np.array([0,2,4,6,8,10])*1e-3
center_flux = 16.6e-6

m = -center_flux / 10

estimated_flux = [m*x + center_flux for x in radial]


delta_r = abs(radial_pos[1] - radial_pos[0])
delta_r = abs(np.mean(np.array(radial_pos[1::]) - np.array(radial_pos[0:-1])))
for count, pos in enumerate(radial_pos):
    if count == 0:
        annulus_area[count] = np.pi * (radial_pos[count] + delta_r / 2.)**2
    else:
        annulus_area[count] = np.pi * ((radial_pos[count] + delta_r / 2.)**2 
            - (radial_pos[count] - delta_r / 2.)**2)

total_spray_area = np.pi*(radial_pos[-1] + delta_r/2)**2
     
print("Estimated liquid flow rate: {:.8f} mL/min".format(np.sum(estimated_flux*annulus_area[0:-1])*60e6))
