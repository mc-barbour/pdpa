# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 13:44:20 2023

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


#%% Load single run data and plot

data_dir = "D:/Barbour/OneDrive - UW/Downstream/Data/"

df_day1 = pd.read_csv(data_dir + '20230112_Downstream/experiment_statistics.csv')
df_day2 = pd.read_csv(data_dir + '20230113_Downstream/experiment_statistics.csv')


#%% caclulate momentum flux ratio

def estimate_flowrate(df):
    return df["Tidal volume (mL)"] / df["Inspiration period (s)"] * 60 / 1000

od_orrifice = 3./32. * 2.54e-2
rho_air = 1.185
rho_water = 1000

M = []
for count in range(len(df_day1)):
    df_tmp = df_day1.loc[count]
    
    if df_tmp['Needle Gauge'] == 22:
        d_needle_liq = 0.016*2.54e-2
        d_needle_air = 0.028*2.54e-2
    elif df_tmp['Needle Gauge'] == 19:
        d_needle_liq = 0.027*2.54e-2
        d_needle_air = 0.042*2.54e-2
        
    A_air = np.pi / 4. * (od_orrifice**2 - d_needle_air**2)
    A_water = np.pi / 4. * d_needle_liq**2

    Q_water = df_tmp['Injection rate (ml/min)']
    V_water = Q_water/1e6/60/A_water

    Q_air = estimate_flowrate(df_tmp)
    V_air = Q_air / 1e3 / 60 / A_air
    
    m = rho_air * V_air**2 / (rho_water * V_water**2)
    M.append(m)
    
df_day1["Momentum ratio"] = M

M = []
for count in range(len(df_day2)):
    df_tmp = df_day2.loc[count]
    
    if df_tmp['Needle Gauge'] == 22:
        d_needle_liq = 0.016*2.54e-2
        d_needle_air = 0.028*2.54e-2
    elif df_tmp['Needle Gauge'] == 19:
        d_needle_liq = 0.027*2.54e-2
        d_needle_air = 0.042*2.54e-2
        
    A_air = np.pi / 4. * (od_orrifice**2 - d_needle_air**2)
    A_water = np.pi / 4. * d_needle_liq**2

    Q_water = df_tmp['Injection rate (ml/min)']
    V_water = Q_water/1e6/60/A_water

    Q_air = estimate_flowrate(df_tmp)
    V_air = Q_air / 1e3 / 60 / A_air
    
    m = rho_air * V_air**2 / (rho_water * V_water**2)
    M.append(m)
    
df_day2["Momentum ratio"] = M




#%% plot velocity, d10, d32, d43 at centerline for different momentum ratios

fig = make_subplots(rows=2, cols=2)

plot_df = df_day2[(df_day2['Radial position (mm)']==0) & (df_day2['Injection rate (ml/min)']==6)]

fig.add_trace(go.Scatter(x=plot_df['Tidal volume (mL)'], y=plot_df['U mean (m/s)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[0],
                         showlegend=False), row=1,col=1)

fig.add_trace(go.Scatter(x=plot_df['Tidal volume (mL)'], y=plot_df['D10 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[0],
                         showlegend=False), row=1,col=2)

fig.add_trace(go.Scatter(x=plot_df['Tidal volume (mL)'], y=plot_df['D32 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[0],
                         showlegend=False), row=2,col=1)

fig.add_trace(go.Scatter(x=plot_df['Tidal volume (mL)'], y=plot_df['D43 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[0],
                         name='19 Gauge Needle'), row=2,col=2)

fig.update_layout(template='plotly_white', font_family='sans_serif', font_size=16)
fig.update_xaxes(title='Tidal Volume Ventilator (mL)')
fig.update_yaxes(title='Mean Velocity, <u> (m/s)', row=1, col=1)
fig.update_yaxes(title='D10 (um)', row=1, col=2)
fig.update_yaxes(title='D32 (um)', row=2, col=1)
fig.update_yaxes(title='D43 (um)', row=2, col=2)

plot_df = plot_df = df_day1[(df_day1['Radial position (mm)']==0) & (df_day1['Injection rate (ml/min)']==6)]
fig.add_trace(go.Scatter(x=plot_df['Tidal volume (mL)'], y=plot_df['U mean (m/s)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[1],
                         showlegend=False), row=1,col=1)

fig.add_trace(go.Scatter(x=plot_df['Tidal volume (mL)'], y=plot_df['D10 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[1],
                         showlegend=False), row=1,col=2)

fig.add_trace(go.Scatter(x=plot_df['Tidal volume (mL)'], y=plot_df['D32 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[1],
                         showlegend=False), row=2,col=1)

fig.add_trace(go.Scatter(x=plot_df['Tidal volume (mL)'], y=plot_df['D43 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[1],
                         name='22 Gauge Needle'), row=2,col=2)

fig.show()

#%% same figure as above but with momentum ratio

fig = make_subplots(rows=2, cols=2)

plot_df = df_day2[(df_day2['Radial position (mm)']==0) & (df_day2['Injection rate (ml/min)']==6)]

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['U mean (m/s)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[0],
                         showlegend=False), row=1,col=1)

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D10 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[0],
                         showlegend=False), row=1,col=2)

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D32 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[0],
                         showlegend=False), row=2,col=1)

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D43 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[0],
                         name='19 Gauge Needle'), row=2,col=2)

fig.update_layout(template='plotly_white', font_family='sans_serif', font_size=16)
fig.update_xaxes(title='Tidal Volume Ventilator (mL)')
fig.update_yaxes(title='Mean Velocity, <u> (m/s)', row=1, col=1)
fig.update_yaxes(title='D10 (um)', row=1, col=2)
fig.update_yaxes(title='D32 (um)', row=2, col=1)
fig.update_yaxes(title='D43 (um)', row=2, col=2)

plot_df = plot_df = df_day1[(df_day1['Radial position (mm)']==0) & (df_day1['Injection rate (ml/min)']==6)]
fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['U mean (m/s)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[1],
                         showlegend=False), row=1,col=1)

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D10 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[1],
                         showlegend=False), row=1,col=2)

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D32 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[1],
                         showlegend=False), row=2,col=1)

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D43 scaled (um)'],
                         mode='markers',
                         marker_size=16,
                         marker_symbol='star-diamond',
                         marker_color=colors[1],
                         name='22 Gauge Needle'), row=2,col=2)

fig.show()

#%% plot velocity, d10, d32, d43 at centerline for different momentum ratios - changing injection flowrate

fig = make_subplots(rows=2, cols=2)

plot_df = df_day2[(df_day2['Radial position (mm)']==0) & (df_day2['Tidal volume (mL)']==850)]

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['U mean (m/s)'],
                          mode='markers',
                          marker_size=16,
                          marker_symbol='star-diamond',
                          marker_color=colors[0],
                          showlegend=False), row=1,col=1)

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D10 scaled (um)'],
                          mode='markers',
                          marker_size=16,
                          marker_symbol='star-diamond',
                          marker_color=colors[0],
                          showlegend=False), row=1,col=2)

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D32 scaled (um)'],
                          mode='markers',
                          marker_size=16,
                          marker_symbol='star-diamond',
                          marker_color=colors[0],
                          showlegend=False), row=2,col=1)

fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D43 scaled (um)'],
                          mode='markers',
                          marker_size=16,
                          marker_symbol='star-diamond',
                          marker_color=colors[0],
                          name='19 Gauge Needle'), row=2,col=2)

fig.update_layout(template='plotly_white', font_family='sans_serif', font_size=16)
fig.update_xaxes(title='Tidal Volume Ventilator (mL)')
fig.update_yaxes(title='Mean Velocity, <u> (m/s)', row=1, col=1)
fig.update_yaxes(title='D10 (um)', row=1, col=2)
fig.update_yaxes(title='D32 (um)', row=2, col=1)
fig.update_yaxes(title='D43 (um)', row=2, col=2)

plot_df = plot_df = df_day1[(df_day1['Radial position (mm)']==0) & (df_day1['Tidal volume (mL)']==750)]
# fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['U mean (m/s)'],
#                          mode='markers',
#                          marker_size=16,
#                          marker_symbol='star-diamond',
#                          marker_color=colors[1],
#                          showlegend=False), row=1,col=1)

# fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D10 scaled (um)'],
#                          mode='markers',
#                          marker_size=16,
#                          marker_symbol='star-diamond',
#                          marker_color=colors[1],
#                          showlegend=False), row=1,col=2)

# fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D32 scaled (um)'],
#                          mode='markers',
#                          marker_size=16,
#                          marker_symbol='star-diamond',
#                          marker_color=colors[1],
#                          showlegend=False), row=2,col=1)

# fig.add_trace(go.Scatter(x=plot_df['Momentum ratio'], y=plot_df['D43 scaled (um)'],
#                          mode='markers',
#                          marker_size=16,
#                          marker_symbol='star-diamond',
#                          marker_color=colors[1],
#                          name='22 Gauge Needle'), row=2,col=2)

fig.show()

#%%Plot by injection flowrate

fig = make_subplots(rows=2, cols=2)

plot_df = df_day2[(df_day2['Radial position (mm)']==0) & (df_day2['Tidal volume (mL)']==850)]

fig.add_trace(go.Scatter(x=plot_df['Injection rate (ml/min)'], y=plot_df['U mean (m/s)'],
                          mode='markers',
                          marker_size=16,
                          marker_symbol='star-diamond',
                          marker_color=colors[0],
                          showlegend=False), row=1,col=1)

fig.add_trace(go.Scatter(x=plot_df['Injection rate (ml/min)'], y=plot_df['D10 scaled (um)'],
                          mode='markers',
                          marker_size=16,
                          marker_symbol='star-diamond',
                          marker_color=colors[0],
                          showlegend=False), row=1,col=2)

fig.add_trace(go.Scatter(x=plot_df['Injection rate (ml/min)'], y=plot_df['D32 scaled (um)'],
                          mode='markers',
                          marker_size=16,
                          marker_symbol='star-diamond',
                          marker_color=colors[0],
                          showlegend=False), row=2,col=1)

fig.add_trace(go.Scatter(x=plot_df['Injection rate (ml/min)'], y=plot_df['D43 scaled (um)'],
                          mode='markers',
                          marker_size=16,
                          marker_symbol='star-diamond',
                          marker_color=colors[0],
                          name='19 Gauge Needle'), row=2,col=2)

fig.update_layout(template='plotly_white', font_family='sans_serif', font_size=16)
fig.update_xaxes(title='Injection Flowrate (mL/min)')
fig.update_yaxes(title='Mean Velocity, <u> (m/s)', row=1, col=1)
fig.update_yaxes(title='D10 (um)', row=1, col=2)
fig.update_yaxes(title='D32 (um)', row=2, col=1)
fig.update_yaxes(title='D43 (um)', row=2, col=2)

    
#%% Calculate full spray metrics

data_dir = '/Users/mbarbour/OneDrive - UW/Downstream/Data/20230103_Downstream/'
data_dir = 'D:/Barbour/OneDrive - UW/Downstream/Data/20230113_Downstream/'

experiment_names = get_experiment_names(data_dir)

experiments_0113 = []

for name in experiment_names:
    experiments_0113.append(load_processed_data(name, data_dir))

data_dir = 'D:/Barbour/OneDrive - UW/Downstream/Data/20230112_Downstream/'

experiment_names = get_experiment_names(data_dir)

experiments_0112 = []

for name in experiment_names:
    experiments_0112.append(load_processed_data(name, data_dir))

#%% execute

params = {"max_diameter": 100,
          "slot_area": 150e-6 / np.sin(150. * np.pi/180) * (750/250),
          "median_particle_diameter": [10,20],
          "number of diameter bins": 51}



pdf_850, VF_850, VFD_850, VFD_full_850, bin_volume_850, total_spray_area_850, bins, diam_dict_850 = process_full_spray(experiments_0113[0:7], params, use_scaled_diameter=True)
pdf_750, VF_750, VFD_750, VFD_full_750, bin_volume_750, total_spray_area_750, bins, diam_dict_750 = process_full_spray(experiments_0112[6:12], params, use_scaled_diameter=True)

#%% generate figures


fig = go.Figure()
fig.add_trace(go.Scatter(x=bins, y=np.cumsum(VFD_850) / np.sum(VFD_850),
                         marker_color=colors[0],
                         name='TV=850mL, QL=6ml/min,<br> 19 G Needle, M=230 '),
              )
fig.add_trace(go.Scatter(x=bins, y=np.cumsum(VFD_750) / np.sum(VFD_750),
                         marker_color=colors[1],
                         name='TV=750mL, QL=6ml/min, <br> 22 G Needle, M=17'),)
fig.update_layout(template='plotly_white', font_family='sans_serif', font_size=16)
fig.update_yaxes(title='CDF(VFDF)')
fig.update_xaxes(title='Particle Diameter (um)')
fig.show()


fig = go.Figure()
fig.add_trace(go.Scatter(x=bins, y=pdf_850,
              marker_color=colors[0],
              name='TV=850mL, QL=6ml/min,<br> 19 G Needle, M=230'))
fig.add_trace(go.Scatter(x=bins, y=pdf_750,
              marker_color=colors[1],
              name='TV=750mL, QL=6ml/min, <br> 22 G Needle, M=17'))
fig.update_yaxes(title='PDF(d)', type='log')
fig.update_layout(template='plotly_white', font_family='sans_serif', font_size=16)
fig.update_xaxes(title='Particle Diameter (um)')


fig.show()

#%% process one of Peter's runs

data_dir = '/Users/mbarbour/OneDrive - UW/Downstream/Data/20230103_Downstream/'
data_dir = 'D:/Barbour/OneDrive - UW/Downstream/Data/peter_data/180SLPM_6mlmin/'

experiment_names = get_experiment_names(data_dir)

experiments= []

for name in experiment_names:
    experiments.append(load_processed_data(name, data_dir))

#%% scale peter's data

for exp in experiments:
    exp.scale_diameter()

stat_dict = single_measurement_statistics(experiments[0], use_scaled_diameters=True)


#%% 240 slpm

params = {"max_diameter": 200,
          "slot_area": 150e-6 / np.sin(150. * np.pi/180) * (750/250),
          "median_particle_diameter": [10,20],
          "number of diameter bins": 51}

pdf_peter, VF_peter, VFD_peter, VFD_full_peter, bin_volume, total_spray_area, bins_peter, diam_dict_peer = process_full_spray(experiments, params, use_scaled_diameter=True)
#%% 180 slpm

params = {"max_diameter": 150,
          "slot_area": 150e-6 / np.sin(150. * np.pi/180) * (750/250),
          "median_particle_diameter": [10,20],
          "number of diameter bins": 51}

pdf_peter_180, VF_peter_180, VFD_peter_180, VFD_full_peter_180, bin_volume_180, total_spray_area, bins_peter_180, diam_dict_peter_180 = process_full_spray(experiments, params, use_scaled_diameter=True)



#%% plot

fig = go.Figure()
fig.add_trace(go.Scatter(x=bins, y=np.cumsum(VFD_850) / np.sum(VFD_850),
                         marker_color=colors[0],
                         mode='lines+markers',
                         name='TV=850mL, QL=6ml/min,<br> 19 G Needle, M=230 '),
              )

fig.add_trace(go.Scatter(x=bins, y=np.cumsum(VFD_750) / np.sum(VFD_750),
                         marker_color=colors[1],
                         mode='lines+markers',
                         name='TV=750mL, QL=6ml/min, <br> 22 G Needle, M=17'),)

fig.add_trace(go.Scatter(x=bins_peter, y=np.cumsum(VFD_peter) / np.sum(VFD_peter),
                         marker_color=colors[2],
                         mode='lines+markers',
                         name='Comp Air: Q_air=240SLPM, <br>W_water=6ml/min, 19 G Needle'),
              )

# fig.add_trace(go.Scatter(x=bins_peter_180, y=np.cumsum(VFD_peter_180) / np.sum(VFD_peter_180),
#                          marker_color=colors[3],
#                          mode='lines+markers',
#                          name='Comp Air: Q_air=180SLPM,<br> W_water=6ml/min, 19 G Needle'),
#               )



fig.update_layout(template='plotly_white', font_family='sans_serif', font_size=16)
fig.update_yaxes(title='CDF(VFDF)')
fig.update_xaxes(title='Particle Diameter (um)')
fig.show()


fig = go.Figure()
fig.add_trace(go.Scatter(x=bins, y=pdf_850,
              marker_color=colors[0],
              mode='lines+markers',
              name='TV=850mL, QL=6ml/min,<br> 19 G Needle')

              )
fig.add_trace(go.Scatter(x=bins, y=pdf_750,
              marker_color=colors[1],
              mode='lines+markers',
              name='TV=750mL, QL=6ml/min, <br> 22 G Needle'))

fig.add_trace(go.Scatter(x=bins_peter, y=pdf_peter,
                         marker_color=colors[2],
                         mode='lines+markers',
                         name='Comp Air: Q_air=240SLPM, <br>W_water=6ml/min, 19 G Needle'),
              )

# fig.add_trace(go.Scatter(x=bins_peter_180, y=pdf_peter_180,
#                          marker_color=colors[3],
#                          mode='lines+markers',
#                          name='Comp Air: Q_air=180SLPM,<br> W_water=6ml/min, 19 G Needle'),)

fig.update_yaxes(title='PDF(d)', type='log')
fig.update_layout(template='plotly_white', font_family='sans_serif', font_size=16,  legend=dict(yanchor="top",
    y=0.99,
    xanchor="right",
    x=0.99))
fig.update_xaxes(title='Particle Diameter (um)', range=(0,120))


fig.show()




