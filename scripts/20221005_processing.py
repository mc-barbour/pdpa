#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Created on Mon Oct 10 10:11:30 2022

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

data_dir = "/Users/mbarbour/OneDrive - UW/Downstream/Data/Neubulizer_v2_testing_20221005/"

experiment_names = sorted(get_experiment_names(data_dir))

experiments = []

for name in experiment_names:

    experiments.append(load_experiment(name, data_dir))
    
#%% scale and filter diameters

for exp in experiments:
    exp.scale_diameter()
    exp.filter_diameters(filter_curves=True)
    
                       
#%% plot droplet diameter vs droplet velocity

nrows=3;ncols=4
fig = make_subplots(rows=nrows, cols=ncols)

for count, exp in enumerate(experiments[0:12]):
    
    if exp.measurement_distance == 30:
        col=int(count/nrows)+1
        print(col)
        fig.add_trace(go.Scatter(x=exp.diameter[::2], y=exp.velocity_chan1[::2], 
                                 mode='markers', name=exp._name,
                                 marker_color=colors[0]), row=1, col=col)
    

        
    elif exp.measurement_distance == 60:
        col=int(count/nrows)+1
        fig.add_trace(go.Scatter(x=exp.diameter[::2], y=exp.velocity_chan1[::2], mode='markers', 
                                 name=exp._name, marker_color=colors[1]), row=2, col=int(count/nrows)+1)
        # fig.update_yaxes(title='Velocity (m/s)', range=[0,50], row=2, col=col)

        
    elif exp.measurement_distance == 90:
        col=int(count/nrows)+1
        fig.add_trace(go.Scatter(x=exp.diameter[::2], y=exp.velocity_chan1[::2], 
                                 mode='markers', name=exp._name,
                                 marker_color=colors[2]), row=3, col=int(count/nrows)+1)
        # fig.update_xaxes(title='Diameter (um), TV = {:g}'.foramt(tidal_volumes), range=[0,800], row=3, col=col)
        # fig.update_yaxes(title='Velocity (m/s)', range=[0,50], row=3, col=col)
fig.update_yaxes(range=[0,50])
fig.update_xaxes(range=[0,100])

fig.update_yaxes(title="Velocity (m/s) <br> 30mm downstream", row=1,col=1)
fig.update_yaxes(title="Velocity (m/s) <br> 60mm downstream", row=2,col=1)
fig.update_yaxes(title="Velocity (m/s) <br> 90mm downstream", row=3,col=1)

fig.update_xaxes(title="Diameter (um), TV=450mL", row=3,col=1)
fig.update_xaxes(title="Diameter (um), TV=550mL", row=3,col=2)
fig.update_xaxes(title="Diameter (um), TV=650mL", row=3,col=3)
fig.update_xaxes(title="Diameter (um), TV=750mL", row=3,col=4)

fig.show()

                       
#%% plot histograms of droplet diameter

nrows=3;ncols=4
fig = make_subplots(rows=nrows, cols=ncols)

for count, exp in enumerate(experiments[0:12]):
    
    if exp.measurement_distance == 30:
        col=int(count/nrows)+1
        print(col)
        fig.add_trace(go.Histogram(x=exp.diameter,
                                   xbins=dict(size=2), autobinx=False,
                                   name=exp._name,
                                   marker_color=colors[0]), row=1, col=col)
    
    elif exp.measurement_distance == 60:
        col=int(count/nrows)+1
        fig.add_trace(go.Histogram(x=exp.diameter,
                                   xbins=dict(size=2), autobinx=False,
                                   name=exp._name, marker_color=colors[1]), row=2, col=int(count/nrows)+1)
        # fig.update_yaxes(title='Velocity (m/s)', range=[0,50], row=2, col=col)

        
    elif exp.measurement_distance == 90:
        col=int(count/nrows)+1
        fig.add_trace(go.Histogram(x=exp.diameter,
                                 name=exp._name,
                                 xbins=dict(size=2), autobinx=False,
                                 marker_color=colors[2]), row=3, col=int(count/nrows)+1)
        # fig.update_xaxes(title='Diameter (um), TV = {:g}'.foramt(tidal_volumes), range=[0,800], row=3, col=col)
        # fig.update_yaxes(title='Velocity (m/s)', range=[0,50], row=3, col=col)
# fig.update_yaxes(range=[0,50])
fig.add_vline(x=10, line_color='grey', line_dash='dash')
fig.update_xaxes(range=[0,100])

fig.update_yaxes(title="Count <br> 30mm downstream", row=1,col=1)
fig.update_yaxes(title="Count <br> 60mm downstream", row=2,col=1)
fig.update_yaxes(title="Count <br> 90mm downstream", row=3,col=1)

fig.update_xaxes(title="Diameter (um), TV=450mL", row=3,col=1)
fig.update_xaxes(title="Diameter (um), TV=550mL", row=3,col=2)
fig.update_xaxes(title="Diameter (um), TV=650mL", row=3,col=3)
fig.update_xaxes(title="Diameter (um), TV=750mL", row=3,col=4)

fig.show()



#%% plot average particle size vs flowrate 

fig = go.Figure()

for exp in experiments:
    
    avg_diameter = [np.nanmean(exp.diameter)]
    print(avg_diameter)
    
    if exp.measurement_distance == 30:
        
        fig.add_trace(go.Scatter(y=avg_diameter, x=[exp.estimated_flowrate], mode='markers',
                                 marker_symbol = 'circle',
                                 marker_size=20,
                                 marker_color=colors[0], showlegend=False))
    
    elif exp.measurement_distance == 60:
        
        fig.add_trace(go.Scatter(y=avg_diameter, x=[exp.estimated_flowrate], mode='markers',
                                 marker_symbol = 'diamond',
                                 marker_size=20,
                                 marker_color=colors[1], showlegend=False))
    
    elif exp.measurement_distance == 90:
        
        fig.add_trace(go.Scatter(y=avg_diameter, x=[exp.estimated_flowrate], mode='markers',
                                 marker_symbol = 'star-triangle-up',
                                 marker_size=20,
                                 marker_color=colors[2], showlegend=False))


fig.add_trace(go.Scatter(x=[0], y=[0], mode='markers',
                         marker_symbol = 'circle',
                         marker_size=20,
                         marker_color=colors[0],
                         name = 'X=30mm'))
fig.add_trace(go.Scatter(x=[0], y=[0], mode='markers',
                         marker_symbol = 'diamond',
                         marker_size=20,
                         marker_color=colors[1],
                         name = 'X=60mm'))
fig.add_trace(go.Scatter(x=[0], y=[0], mode='markers',
                         marker_symbol = 'star-triangle-up',
                         marker_size=20,
                         marker_color=colors[2],
                         name = 'X=90mm'))


fig.update_xaxes(title='Ventilator Flowrate (LPM)', range=[12,24])
fig.update_yaxes(title='Average Diameter (um)', range=[13,25])

fig.show()


#%% plot average particle size 

fig = go.Figure()

for exp in experiments:
    
    avg_diameter = [np.nanmean(exp.diameter)]
    d20 = (np.nansum(exp.diameter**2) / np.nansum(exp.diameter**0))**(1 / 2)
    d30 = (np.nansum(exp.diameter**3) / np.nansum(exp.diameter**0))**(1 / 3)
    d32 = (np.nansum(exp.diameter**3) / np.nansum(exp.diameter**2))**(1 / 1)
    d43 = (np.nansum(exp.diameter**4) / np.nansum(exp.diameter**3))**(1 / 1)
    
    
    print(avg_diameter)
    
    if exp.tidal_volume == 450:
        
        fig.add_trace(go.Scatter(y=avg_diameter, x=[exp.measurement_distance], mode='markers',
                                 marker_symbol = 'circle',
                                 marker_size=20,
                                 marker_color=colors[0], showlegend=False))
    elif exp.tidal_volume == 550:
       
       fig.add_trace(go.Scatter(y=avg_diameter, x=[exp.measurement_distance], mode='markers',
                                marker_symbol = 'hexagon',
                                marker_size=20,
                                marker_color=colors[1], showlegend=False))
    
    elif exp.tidal_volume == 650:
        
        fig.add_trace(go.Scatter(y=avg_diameter, x=[exp.measurement_distance], mode='markers',
                                 marker_symbol = 'star-triangle-up',
                                 marker_size=20,
                                 marker_color=colors[2], showlegend=False))
    elif exp.tidal_volume == 750:
        
        fig.add_trace(go.Scatter(y=avg_diameter, x=[exp.measurement_distance], mode='markers',
                                 marker_symbol = 'star-triangle-down',
                                 marker_size=20,
                                 marker_color=colors[3], showlegend=False))



fig.add_trace(go.Scatter(x=[0], y=[0], mode='markers',
                         marker_symbol = 'circle',
                         marker_size=20,
                         marker_color=colors[0],
                         name = 'TV=450'))
fig.add_trace(go.Scatter(x=[0], y=[0], mode='markers',
                         marker_symbol = 'hexagon',
                         marker_size=20,
                         marker_color=colors[1],
                         name = 'TV=550'))
fig.add_trace(go.Scatter(x=[0], y=[0], mode='markers',
                         marker_symbol = 'star-triangle-up',
                         marker_size=20,
                         marker_color=colors[2],
                         name = 'TV=650'))
fig.add_trace(go.Scatter(x=[0], y=[0], mode='markers',
                         marker_symbol = 'star-triangle-up',
                         marker_size=20,
                         marker_color=colors[3],
                         name = 'TV=750'))


fig.update_xaxes(title='Ventilator Flowrate (LPM)', range=[10,100])
fig.update_yaxes(title='Average Diameter (um)', range=[13,25])

fig.show()



#%% Test for calculating distribution function

nbins = 5000

bins = np.linspace(0,150,1000)

fig = make_subplots(rows=2, cols=1)
for exp in experiments:
    
    
    diameters = exp.diameter[~np.isnan(exp.diameter)]
    
    
    # pdf_V, bins_V = np.histogram(exp.volume, bins = nbins, range=(0,50), density=True)
    pdf_D, bins_D = np.histogram(exp.diameter, bins=nbin, range=(0,100), density=True)
    pdf_u, bins_u = np.histogram(exp.velocity_chan1, bins=nbin, density=True)
    # bins_centers_V = 0.5 * (bins_V[1:] + bins_V[:-1])
    bins_centers_D = 0.5 * (bins_D[1:] + bins_D[:-1])
    bins_centers_u = 0.5 * (bins_u[1:] + bins_u[:-1])
    
    


    fig.add_trace(go.Scatter(x = bins_centers_D, y=pdf_D, name=exp._name), row=1, col=1)
    fig.add_trace(go.Scatter(x = bins_centers_V, y=pdf_V, name=exp._name), row=2, col=1)

    # fig.add_trace(go.Scatter(x = 2. * (3/4*np.pi*bins_centers_V)**(1/3),
    #                          y = pdf_V, name = exp._name))

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



#%%

nbins = 75
diameter_threshold = 150

bins = np.linspace(0, diameter_threshold, nbins)

fig = make_subplots(rows=3, cols=1)

plot_count = 0
for count, exp in enumerate(experiments):
    
 
    if exp.measurement_distance == 90:
        diameters = exp.diameter[~np.isnan(exp.diameter)]
        diameters = diameters[np.where(diameters < diameter_threshold)[0]]
        
        bincount = np.bincount(np.digitize(diameters,bins[0:-1]), minlength=nbins)
        
        bin_volume = 4. / 3. * np.pi * (bins / 2.) **3 * bincount # volume in each bin
        
        bin_percent_count = bincount * (1 / np.sum(bincount)) * 100
        bin_percent_volume = bin_volume * (1 / np.sum(bin_volume)) * 100
        
        bin_count_cum = np.cumsum(bin_percent_count, dtype=float)
        bin_volume_cum = np.cumsum(bin_percent_volume, dtype=float)
        
        v90 = intersection(bin_volume_cum, 90, bins)
        print(exp._name, v90)
        
        fig.add_trace(go.Scatter(x=bins, y=bincount, name=exp._name, marker_color=colors[plot_count]), row=1, col=1)
        
        fig.add_trace(go.Scatter(x=bins, y=bin_percent_volume, name=exp._name, marker_color=colors[plot_count]), row=2, col=1)
    
        fig.add_trace(go.Scatter(x=bins, y=bin_volume_cum, name=exp._name, marker_color=colors[plot_count]), row=3, col=1)
        fig.add_vline(v90, line_color=colors[plot_count], row=3, col=1)
        plot_count = plot_count + 1
        
        

fig.update_yaxes(title='bin counts', row=1, col=1)
fig.update_xaxes(title='diameter', row=1, col=1)

fig.update_yaxes(title='volume percentage', row=2, col=1)
fig.update_xaxes(title='diameter', row=2, col=1)
    
fig.update_yaxes(title='volume cumsum', row=3, col=1)
fig.update_xaxes(title='diameter', row=3, col=1)


fig.show()

#%% compute v90 for all runs

nbins = 500
diameter_threshold = 150

bins = np.linspace(0,diameter_threshold,nbins)

fig = make_subplots(rows=3, cols=1)

plot_count = 0

v90 = []
measurement_distance = []
tidal_volume = []

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
    
        

df = pd.DataFrame(data=np.array([tidal_volume, measurement_distance, v90]).T,
                  columns=['Tidal Volume', 'Measurement Distance', 'v90'])


fig = go.Figure()

fig.add_trace(go.Scatter(x=df[df['Measurement Distance']==90]['Tidal Volume'],
                         y=df[df['Measurement Distance']==90].v90,
                         name = '90mm',
                         mode = 'markers', 
                         marker_color=colors[0],
                         marker_size=20,
                         marker_symbol='circle')
                          
              )
fig.add_trace(go.Scatter(x=df[df['Measurement Distance']==60]['Tidal Volume'],
                         y=df[df['Measurement Distance']==60].v90,
                         name = '60mm',
                         mode = 'markers', 
                         marker_color=colors[1],
                         marker_size=20,
                         marker_symbol='hexagon')
                          
              )
fig.add_trace(go.Scatter(x=df[df['Measurement Distance']==30]['Tidal Volume'],
                         y=df[df['Measurement Distance']==30].v90,
                         name = '30mm',
                         mode = 'markers', 
                         marker_color=colors[2],
                         marker_size=20,
                         marker_symbol='star-triangle-up')
                          
              )
fig.update_yaxes(title='D90 (um)')
fig.update_xaxes(title='Target tidal volume (mL)')


fig.show()

#%% diameter difference vs diameter


for exp in experiments[0:2]:
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=exp.diameter, y=exp.diameter_diff,
                             mode='markers'))
    fig.update_yaxes(title='Diameter difference (um)')
    fig.update_xaxes(title='Diameter (um)')

    fig.show()

