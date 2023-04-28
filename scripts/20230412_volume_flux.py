# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 12:25:01 2023

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


#%% Load processed data
experiment_names = get_experiment_names(data_dir)

experiments = []

for name in experiment_names:
    experiments.append(load_processed_data(name, data_dir))
    



    
#%%

params = {"max_diameter": 100,
          "slot_width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "median_particle_diameter": [10,20],
          "number of diameter bins": 51}

f, average_beam_length, bin_beam_length, bin_travel_time = calc_beam_length_function(experiments[0:4], params, use_scaled_diameter=False, visualize=True)



#%% let's revisit the beam length

path_length_cutoff = 0.95

use_scaled_diameter = False
measurement_width = params['slot_width']
max_diameter = params['max_diameter']
median_particle_diameter = params["median_particle_diameter"]
num_dbin = params["number of diameter bins"]

measurement_area = measurement_width*1e-3 # initial guess of the measurement area
### Define spray radius and annulus area
radial_pos, annulus_area, total_spray_area = compute_annulus_area(experiments)

bins = np.linspace(0, max_diameter, num_dbin)
bin_centers = bins + (bins[1] - bins[0]) / 2.
bin_centers = bin_centers[0:-1]

# create data structures
bin_particle_count = np.zeros(len(bins))
beam_length_sum = np.zeros(len(bins))
bin_travel_time = np.zeros(len(bins))
average_diameter = 0
total_particle_count = 0
average_beam_length = 0

# calculate the diameter depentend beam length, average beam length, and particle size
for count, exp in enumerate(experiments[0:4]):
    valid_int = np.where(exp.valid == 1)[0]
    
    if use_scaled_diameter:
        d = exp.scaled_diameter[valid_int]
    else:
        d = exp.diameter[valid_int]
        
    bin_count, bin_edge = np.histogram(d, bins)
    digitized_diameter = np.digitize(d, bins)
    
    lx = abs(exp.gate_time[valid_int] * 1e-6 * exp.velocity_chan1[valid_int]) #beam length for each particle
    tx = exp.gate_time[valid_int] * 1e-6
    
    for j in range(num_dbin - 1):
        bin_ind = np.where(digitized_diameter==j+1)[0]
        if len(bin_ind) > 0:
            n_particles_bin = len(bin_ind)
        else:
            n_particles_bin = 0
        
        bin_particle_count[j] = bin_particle_count[j] + n_particles_bin * annulus_area[count] / measurement_area
        beam_length_sum[j] = beam_length_sum[j] + np.nansum(lx[bin_ind]) * annulus_area[count] / measurement_area
        bin_travel_time[j] = bin_travel_time[j] + np.nansum(tx[bin_ind]) * annulus_area[count] / measurement_area
        
        average_diameter = average_diameter + np.nansum(d[bin_ind]) * annulus_area[count] / measurement_area
        total_particle_count = total_particle_count + n_particles_bin * annulus_area[count] / measurement_area
        average_beam_length = average_beam_length + np.nansum(lx[bin_ind]) * annulus_area[count] / measurement_area

non_zero_bins = np.where(bin_particle_count != 0)[0]
bin_beam_length = beam_length_sum / bin_particle_count
bin_travel_time = bin_travel_time / bin_particle_count
average_diameter = average_diameter / total_particle_count
average_beam_length = average_beam_length / total_particle_count

#%% let's really revisit the beam length 

num_dbin = 51

path_length_cutoff = 99
bins = np.linspace(0, max_diameter, num_dbin)

diam_dependent_beamlength = np.empty(num_dbin, dtype=object)

fig=go.Figure()
for count, exp in enumerate(experiments[0:4]):
    # diam_dependent_beamlength = np.empty(num_dbin, dtype=object)

    valid_int = np.where(exp.valid == 1)[0]
    
    d = exp.diameter[valid_int]
        
    bin_count, bin_edge = np.histogram(d, bins)
    digitized_diameter = np.digitize(d, bins)
    
    lx = abs(exp.gate_time[valid_int] * 1e-6 * exp.velocity_chan1[valid_int]) #beam length for each particle
    tx = exp.gate_time[valid_int] * 1e-6
    
    for j in range(num_dbin - 1):
        bin_ind = np.where(digitized_diameter==j+1)[0]
        if len(bin_ind) > 0:
            n_particles_bin = len(bin_ind)
        else:
            n_particles_bin = 0
    
        if count == 0 :
            diam_dependent_beamlength[j] = np.array(lx[bin_ind])
        else:
            diam_dependent_beamlength[j] = np.append(diam_dependent_beamlength[j], lx[bin_ind])

path_length = np.zeros(len(bins))

fig=go.Figure()
for j in range(num_dbin - 1):
    if len(diam_dependent_beamlength[j])==0:
        continue
    path_length[j] = np.percentile(diam_dependent_beamlength[j], path_length_cutoff)

fig.add_trace(go.Scatter(x=bins, y=path_length, name = count, mode='markers'))
fig.update_yaxes(title='Path length (90% percentile) (m)')
fig.update_xaxes(title='Diameter class size (um)')
fig.show()


#%% plot normalized data

d_max = 50

max_index = np.argmin(abs(d_max - bins))
normalized_path_length = np.array(path_length[0:max_index]) / path_length[max_index]
normalized_diameter = bins[0:max_index] / bins[max_index] 

fig = go.Figure()
fig.add_trace(go.Scatter(x = normalized_diameter, y=normalized_path_length))
fig.update_yaxes(title=r'Normalized Beam Length $L(d_i)_{90\%} / L(d_{max})_{90\%}$')
fig.update_xaxes(title=r'$Normalized Diameter (d_i / d_{max})$')
fig.show()




#%% Can we fit a curve to this

def beam_length_func(x, a, b, c, d):
    # return a + b * x / (1 + c * np.exp(0.005 * x))
    return a + b*x + c*x**2 + d*x**3

popt, pcov = curve_fit(beam_length_func, bins[0:17], path_length[0:17], bounds=(0, [np.inf, np.inf, np.inf, np.inf]))

beam_length_fit = beam_length_func(bins, popt[0], popt[1], popt[2], popt[3])

fig = go.Figure()
fig.add_trace(go.Scatter(x=bins, y=path_length, name = count, mode='markers'))
fig.add_trace(go.Scatter(x=bins, y=beam_length_fit, name = count, mode='markers'))
fig.update_yaxes(title='Path length (90% percentile) (m)')
fig.update_xaxes(title='Diameter class size (um)')
fig.show()

    
#%% plot scatter plots of beam length

fig = go.Figure()

for j in range(num_dbin - 1):
    fig.add_trace(go.Scatter(x=np.ones(len(diam_dependent_beamlength[j]))*j, y=diam_dependent_beamlength[j], mode='markers'))
fig.show()


fig = go.Figure()

for j in range(num_dbin - 1):
    fig.add_trace(go.Box(y=diam_dependent_beamlength[j], marker_color=colors[0]))
fig.update_yaxes(title='path length (m)')
fig.show()


#%% percentile

path_length = []

for j in range(num_dbin - 1):
    if len(diam_dependent_beamlength[j])==0:
        continue
    path_length.append(np.percentile(diam_dependent_beamlength[j], 90))


fig=go.Figure()
fig.add_trace(go.Scatter(y=path_length))
fig.show()


#%% okay, let's calculate the volume flux at each location - let's do it in time - we need to split this up by burst

params = {"max_diameter": 60,   # based off of beam width
          "slot_width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "median_particle_diameter": [10,20],
          "number of diameter bins": 51}

min_burst_sep = 0.5
measurement_width = params['slot_width']
max_diameter = params['max_diameter']
median_particle_diameter = params["median_particle_diameter"]
num_dbin = params["number of diameter bins"]
# average_beam_length = params['average beam length']
# bin_travel_time = params['bin travel time']
# f = params['beam length function']

bins = np.linspace(0, max_diameter, num_dbin)
bin_centers = bins + (bins[1] - bins[0]) / 2.
bin_centers = bin_centers[0:-1]

### Define spray radius and annulus area
radial_pos, annulus_area, total_spray_area = compute_annulus_area(experiments)

 
# compute spray parameters
count = 0
A = average_beam_length * measurement_width #slit * np.sin(receiver_angle * np.pi / 180) * 750. / 250.
V = A * average_beam_length

# VFD = np.zeros((len(annulus_area), len(bin_centers)))
# VF = np.zeros((len(annulus_area), len(bin_centers)))
# NFD = np.zeros(len(bins)-1)

# we 

# for count, exp in enumerate(experiments):

exp = experiments[0]
time = exp.time_chan1[valid_int]
valid_int = np.where(exp.valid == 1)[0]
d = exp.diameter[valid_int]
tx = exp.gate_time[valid_int] * 1e-6

burst_seperator = np.where(np.diff(time) > min_sep)[0]
VFD = np.zeros((len(burst_seperator)-1, len(bin_centers)))

for count in range(len(burst_seperator)-1):
    print(count)
    if count == 0:
        burst_ind = np.arange(0, burst_seperator[count])
    else:
        burst_ind = np.arange(burst_seperator[count]+1,burst_seperator[count+1])
        
    burst_time = max(time[burst_ind] - min(time[burst_ind]))
    print(burst_time)
    
    d_burst = d[burst_ind]
    bin_count, bin_edge = np.histogram(d_burst, bins)
    digitized_diameter = np.digitize(d_burst, bins)
    
    for d_bin in range(num_dbin - 1):
        
        corrected_beam_length = beam_length_func(bin_centers[d_bin], popt[0], popt[1], popt[2], popt[3])
       
        beam_area = measurement_width * corrected_beam_length
        bin_ind = np.where(digitized_diameter==d_bin+1)[0]
        
        # print(corrected_beam_length, d_bin)
        if len(bin_ind) == 0:
            tx_bin_mean = 0
            bin_particle_count = 0
        else:
            tx_bin_mean = bin_travel_time[d_bin]
            bin_particle_count = len(bin_ind)
        
        VFD[count, d_bin] = np.pi / 6 * (bin_particle_count * (bins[d_bin]*1e-6)**3) / beam_area / burst_time
       
 
       # NFD[d_bin] = NFD[d_bin] + bin_particle_count / T / A / (1+f[d_bin])*annulus_area[count] # area might be wrong        
        
        # VF[count, d_bin] = bin_particle_count * (bins[d_bin]*1e-6)**3 * tx_bin_mean / T / V / (1.+f[d_bin])

print(np.mean(np.sum(VFD,axis=1)))
#%% 

fig = go.Figure()
fig.add_trace(go.Scatter(y = np.sum(VFD, axis=1)))
fig.update_yaxes(title="Volume Flux Density (m3/( m2 s))")
fig.update_xaxes(title ="Burst / spray event")
fig.show()

#%% now let's look at the volume flex density across the spray

def calc_spray_time(time, d, min_sep, small_count=3, debug=True):

    bad_points = []
    burst_seperator = np.where(np.diff(time) > min_sep)[0]
    print(burst_seperator)
    for count, point in enumerate(burst_seperator):
        if (time[point] - time[point-small_count]) >= (min_sep / 2):
            bad_points.append(count)
    
    burst_seperator = np.delete(burst_seperator, bad_points)
    print(burst_seperator)
    burst_time = []
    for count in range(len(burst_seperator)):
 
        time_chunk = time[0:burst_seperator[count]]
        reverse_index = np.where(abs(np.diff(np.flip(time_chunk)))>1)[0]
        if (count == 0) and (len(reverse_index)==0):
            burst_ind = np.arange(0, burst_seperator[count])
        else:
            burst_start = len(time_chunk) - reverse_index[0]
            burst_ind = np.arange(burst_start,burst_seperator[count])
   
        dt = max(time[burst_ind] - min(time[burst_ind]))
        burst_time.append(dt)
    

    # check for last burst
    if time[burst_seperator[-1]] < len(time):
            burst_ind = np.arange(burst_seperator[-1]+1, len(time))
            dt = max(time[burst_ind] - min(time[burst_ind]))
            burst_time.append(dt)
    print(burst_time)
    print("Average Droplet Burst: {:.4f}s \n Total spray time: {:.4f}s".format(np.mean(burst_time), np.sum(burst_time)))
    if debug:
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=time, y=d))
        fig.add_trace(go.Scatter(x=time[burst_seperator], y=d[burst_seperator], marker_color='red',mode='markers'))
        
        fig.show()
    
    return np.sum(burst_time)

exp = experiments[2]
valid_int = np.where(exp.valid == 1)[0]
d = exp.diameter[valid_int]
time = exp.time_chan1[valid_int]
spray_time = calc_spray_time(time, d, 0.5)




#%% going to divide by the burst spray time

meas_locations = len(experiments[0:4])
VFD = np.zeros((meas_locations, len(bin_centers)))

for count, exp in enumerate(experiments[0:4]):
    
    valid_int = np.where(exp.valid == 1)[0]
    time = exp.time_chan1[valid_int]
    d = exp.diameter[valid_int]
    
    spray_time = calc_spray_time(time, d, min_burst_sep)
    
    valid_int = np.where(exp.valid == 1)[0]
   
    tx = exp.gate_time[valid_int] * 1e-6 
    
    bin_count, bin_edge = np.histogram(d, bins)
    digitized_diameter = np.digitize(d, bins)
    
    for d_bin in range(num_dbin - 1):
        
        corrected_beam_length = beam_length_func(bin_centers[d_bin], popt[0], popt[1], popt[2], popt[3])
       
        beam_area = measurement_width * corrected_beam_length
        bin_ind = np.where(digitized_diameter==d_bin+1)[0]
        
        # print(corrected_beam_length, d_bin)
        if len(bin_ind) == 0:
            tx_bin_mean = 0
            bin_particle_count = 0
        else:
            tx_bin_mean = bin_travel_time[d_bin]
            bin_particle_count = len(bin_ind)
        
        VFD[count, d_bin] = np.pi / 6 * (bin_particle_count * (bins[d_bin]*1e-6)**3) / beam_area / spray_time


#%% 


fig = go.Figure()
fig.add_trace(go.Scatter(x=np.arange(4)*2, y=np.sum(VFD, axis=1)))
fig.update_xaxes(title='Radial Position')
fig.update_yaxes(title="Volume Flux Density (m3/( m2 s))")
fig.show()




#%%
small_count = 3 
debug=True
burst_seperator = np.where(np.diff(time) > min_sep)[0]

bad_points = []

for count, point in enumerate(burst_seperator):
    if (time[point] - time[point-small_count]) >= (min_sep / 2):
        bad_points.append(count)

burst_seperator = np.delete(burst_seperator, bad_points)

burst_time = []
for count in range(len(burst_seperator)-1):
     
    time_chunk = time[0:burst_seperator[count]]
    reverse_index = np.where(abs(np.diff(np.flip(time_chunk)))>1)[0]
    if (count == 0) and (len(reverse_index)==0):
        burst_ind = np.arange(0, burst_seperator[count])
    else:
        burst_start = len(time_chunk) - reverse_index[0]
        burst_ind = np.arange(burst_start,burst_seperator[count])
   
    dt = max(time[burst_ind] - min(time[burst_ind]))
    burst_time.append(dt)
print(burst_time)
    
print("Average Droplet Burst: {:.4f}s \n Total spray time: {:.4f}s".format(np.mean(burst_time), np.sum(burst_time)))
if debug:
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=time, y=d))
    fig.add_trace(go.Scatter(x=time[burst_seperator], y=d[burst_seperator], marker_color='red',mode='markers'))

    fig.show()


#%% let's create a new function for the beam length
    
def beam_length_func_poly3(x, a, b, c, d):

    return a + b*x + c*x**2 + d*x**3

def beam_length_func_poly2(x, a, b, c):

    return a + b*x + c*x**2 

def beam_length_func_plateau(x, a, b, c):
    return a + b * x / (1 + c * np.exp(0.005 * x))



def calc_diameter_beam_length(experiments, params, visualize=True):
    
    
    max_diameter = params['max diameter']
    num_dbin = params["number of diameter bins"]   
    path_length_cutoff = params["beam length cutoff"]
    function_type = params['function type']
    
    bins = np.linspace(0, max_diameter, num_dbin)
    diam_dependent_beamlength = np.empty(num_dbin, dtype=object)
    

    for count, exp in enumerate(experiments[0:4]):
        # diam_dependent_beamlength = np.empty(num_dbin, dtype=object)
    
        valid_int = np.where(exp.valid == 1)[0]
        
        d = exp.diameter[valid_int]
            
        bin_count, bin_edge = np.histogram(d, bins)
        digitized_diameter = np.digitize(d, bins)
        
        lx = abs(exp.gate_time[valid_int] * 1e-6 * exp.velocity_chan1[valid_int]) #beam length for each particle
        tx = exp.gate_time[valid_int] * 1e-6
        
        for j in range(num_dbin - 1):
            bin_ind = np.where(digitized_diameter==j+1)[0]
            if len(bin_ind) > 0:
                n_particles_bin = len(bin_ind)
            else:
                n_particles_bin = 0
        
            if count == 0 :
                diam_dependent_beamlength[j] = np.array(lx[bin_ind])
            else:
                diam_dependent_beamlength[j] = np.append(diam_dependent_beamlength[j], lx[bin_ind])
    
    
    path_length = np.zeros(len(bins))
    for j in range(num_dbin - 1):
        if len(diam_dependent_beamlength[j])==0:
            continue
        path_length[j] = np.percentile(diam_dependent_beamlength[j], path_length_cutoff)
        
        
    if function_type == 'poly3':
        popt, pcov = curve_fit(beam_length_func_poly3, bins[0:-1], path_length[0:-1], bounds=(0, [np.inf, np.inf, np.inf, np.inf]))
        beam_length_fit = beam_length_func_poly3(bins, popt[0], popt[1], popt[2], popt[3])
    
    elif function_type == 'poly2':
        popt, pcov = curve_fit(beam_length_func_poly2, bins[0:-1], path_length[0:-1], bounds=(0, [np.inf, np.inf, np.inf]))
        beam_length_fit = beam_length_func_poly2(bins, popt[0], popt[1], popt[2])
        
    elif function_type == 'plateau':
        popt, pcov = curve_fit(beam_length_func_plateau, bins[0:-1], path_length[0:-1], bounds=(0, [np.inf, np.inf, np.inf]))
        beam_length_fit = beam_length_func_plateau(bins, popt[0], popt[1], popt[2])

    if visualize:
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=bins, y=path_length, name = count, mode='markers'))
        fig.update_yaxes(title='Path length (90% percentile) (m)')
        fig.update_xaxes(title='Diameter class size (um)')
        fig.show()
        
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=bins, y=path_length, name = count, mode='markers'))
        fig.add_trace(go.Scatter(x=bins, y=beam_length_fit, name = count, mode='markers'))
        fig.update_yaxes(title='Path length (cutoff) (m)')
        fig.update_xaxes(title='Diameter class size (um)')
        fig.show()    
            
    return popt
    
params = {"max diameter": 54,
          "slot width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "number of diameter bins": 51,
          "beam length cutoff": 95,
          "function type": 'plateau'}
 
beam_length_func_params = calc_diameter_beam_length(experiments[0:4], params, visualize=True)


#%%

def compute_beam_length(diam, params, popt):
    
    function_type = params['function type']
    
    if function_type == 'poly3':
        beam_length = beam_length_func_poly3(diam, popt[0], popt[1], popt[2], popt[3])
    
    elif function_type == 'poly2':
        beam_length = beam_length_func_poly2(diam, popt[0], popt[1], popt[2])
        
    elif function_type == 'plateau':
        beam_length = beam_length_func_plateau(diam, popt[0], popt[1], popt[2])
        
    return beam_length

def full_spray_metrics(experiments, beam_length_func_params, params):
    
    max_diameter = params['max diameter']
    num_dbin = params["number of diameter bins"]   
    min_burst_sep = params["minimum burst seperation"]
    spray_mode = params['spray mode']
    
    bins = np.linspace(0, max_diameter, num_dbin)
    bin_centers = bins + (bins[1] - bins[0]) / 2.
    bin_centers = bin_centers[0:-1]
    
    ### Define spray radius and annulus area
    radial_pos, annulus_area, total_spray_area = compute_annulus_area(experiments)
    
    
    meas_locations = len(experiments)
    volume_flux = np.zeros((meas_locations, len(bin_centers))) 
    volume_fraction = np.zeros((meas_locations, len(bin_centers)))
    number_density = np.zeros(len(bin_centers))
    
    for count, exp in enumerate(experiments):
        
        valid_int = np.where(exp.valid == 1)[0]
        time = exp.time_chan1[valid_int]
        d = exp.diameter[valid_int]
        tx = exp.gate_time[valid_int] * 1e-6 
        
        if spray_mode == 'burst':
            spray_time = calc_spray_time(time, d, min_burst_sep)
        elif spray_mode == 'continuous':
            spray_time = max(time)
        
        
        bin_count, bin_edge = np.histogram(d, bins)
        digitized_diameter = np.digitize(d, bins)
        
        for d_bin in range(num_dbin - 1):
            
            diam_class_beam_length = compute_beam_length(bin_centers[d_bin], params, beam_length_func_params)
           
            beam_area = measurement_width * diam_class_beam_length
            bin_ind = np.where(digitized_diameter==d_bin+1)[0]  #### check this!!!
            
            if len(bin_ind) == 0:
                tx_bin_mean = 0
                bin_particle_count = 0
            else:
                tx_bin_mean = bin_travel_time[d_bin]
                bin_particle_count = len(bin_ind)
            
            volume_flux[count, d_bin] = np.pi / 6 * (bin_particle_count * (bins[d_bin]*1e-6)**3) / beam_area / spray_time
            # volume_fraction[count, d_bin] = bin_particle_count * (bins[d_bin]*1e-6)**3 * tx_bin_mean / spray_time / V / (1.+f[d_bin])
     
            number_density[d_bin] = number_density[d_bin] + bin_particle_count / spray_time / beam_area * annulus_area[count] # I'm pressure area just gets integrated out of all quantities here       


    ### Compute Metrics to return
    
    pdfN = number_density / np.trapz(number_density, x=bin_centers)
    bin_volume = number_density * np.pi / 6 * (bins[0:-1]*1e-6)**3
    D10T=np.trapz(bin_centers**1*pdfN, x=bin_centers)**(1./1);
    D30T=np.trapz(bin_centers**3*pdfN, x=bin_centers)**(1./3);
    D32T=np.trapz(bin_centers**3*pdfN, x=bin_centers)/np.trapz(bin_centers**2*pdfN, x=bin_centers);
    D43T=np.trapz(bin_centers**4*pdfN, x=bin_centers)/np.trapz(bin_centers**3*pdfN, x=bin_centers);

    diam_dict = {"D10_full": D10T,
                 "D30_full": D30T,
                 "D32_full": D32T,
                 "D43_full": D43T};


    ### integrate VFD and VF across annulus
        
    # volume_fraction_integrated = np.zeros(len(bins)-1)
    volume_flux_integrated = np.zeros(len(bins)-1)
    for bin_count in range(num_dbin - 1):
        # volume_fraction_integrated[bin_count] = np.sum(VF[:, bin_count] * annulus_area)
        volume_flux_integrated[bin_count] = np.sum(VFD[:, bin_count] * annulus_area)

    return pdfN, volume_flux_integrated, volume_flux, bin_volume, total_spray_area, bins, diam_dict



params = {"max diameter": 54,
          "slot width": 150e-6 / np.sin(120. * np.pi/180) * (750/250),
          "number of diameter bins": 51,
          "beam length cutoff": 95,
          "function type": 'plateau', 
          "spray mode": 'burst',
          "minimum burst seperation": 1.0}

pdfN, volume_flux_spray, volume_flux_position, bin_volume, total_spray_area, bins, diam_dict = full_spray_metrics(experiments[0:4], beam_length_func_params, params)

fig = go.Figure()
fig.add_trace(go.Scatter(y=np.sum(volume_flux_position, axis=1)))
fig.show()

#%% okay, lets clear everything and trying running once more
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

fig = go.Figure()
fig.add_trace(go.Scatter(y=np.sum(volume_flux_position, axis=1)))
fig.show()