# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 15:05:35 2023

@author: barbourm
"""

import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
pio.renderers.default = "browser"


def compute_annulus_area(experiments): 
    """
    Calculate annnuls areas for each radial measurement position
    Assumes even spacing for each radial measurement
    """
    radial_pos = []
    
    for exp in experiments:
        radial_pos.append(exp.radial_position * 1e-3)
        
    num_rad_pos = len(radial_pos)
    annulus_area = np.zeros(num_rad_pos)
    
    
    delta_r = abs(radial_pos[1] - radial_pos[0])
    delta_r = abs(np.mean(np.array(radial_pos[1::]) - np.array(radial_pos[0:-1])))
    for count, pos in enumerate(radial_pos):
        if count == 0:
            annulus_area[count] = np.pi * (radial_pos[count] + delta_r / 2.)**2
        else:
            annulus_area[count] = np.pi * ((radial_pos[count] + delta_r / 2.)**2 
                - (radial_pos[count] - delta_r / 2.)**2)
    
    total_spray_area = np.pi*(radial_pos[-1] + delta_r/2)**2
         
    print(np.sum(annulus_area), total_spray_area)   

    return radial_pos, annulus_area, total_spray_area


def intersection(values, cutoff, BINS):

    counter = np.argmax(values >= cutoff) # find cutoff point
    point2 = np.array([BINS[counter], values[counter]]) # cutoff point
    point1 = np.array([BINS[counter-1], values[counter-1]]) # prior point
    slope = (point2[1] - point1[1])/(point2[0] - point1[0])
    intercept = point2[1] - slope*point2[0]
    dist = (cutoff - intercept) * (1/slope) # bin value exactly at cutoff - interpolated
    return dist


def diameter_volume_single_meas(exp, nbins=100, volume_percent=90):
    
    # Calculate the dvXX where volume percent is XX (3g: Dv90)
    
    bins = np.linspace(0, exp.maxdiameter, nbins)
    valid_int = np.where(exp.valid == 1)[0]
    diameters = exp.scaled_diameter[valid_int]
    
    bincount = np.bincount(np.digitize(diameters,bins[0:-1]), minlength=nbins)
    
    bin_volume = 4. / 3. * np.pi * (bins / 2.) **3 * bincount # volume in each bin
    
    bin_percent_count = bincount * (1 / np.sum(bincount)) * 100
    bin_percent_volume = bin_volume * (1 / np.sum(bin_volume)) * 100
    
    bin_count_cum = np.cumsum(bin_percent_count, dtype=float)
    bin_volume_cum = np.cumsum(bin_percent_volume, dtype=float)
    
    diam_volume = intersection(bin_volume_cum, volume_percent, bins)

    return diam_volume

def single_measurement_statistics(exp, use_scaled_diameters=True):
    """
    Calcualte mean and RMS velocity and diameter statistics for a single PDPA data collection run.
    If Scaled_diameter = True
    also calculates diameter statistics using diameter scaled by surface tension ratio
    
    """
    valid_ind = np.where(exp.valid == 1)[0]
    
    if len(valid_ind) == 0:
        print("Warning! All particles are considered valid. Have you performed filtering!?")
    
    characteristics_dict = {
        "U mean (m/s)": np.nanmean(exp.velocity_chan1[valid_ind]),
        "U rms (m/s)" : np.nanstd(exp.velocity_chan1[valid_ind]),
        'D10 (um)'    : np.nanmean(exp.diameter[valid_ind]),
        'D20 (um)'    : np.nansum(exp.diameter[valid_ind]**2 / np.nansum(exp.diameter[valid_ind]**0))**(1./2),
        'D30 (um)'    : np.nansum(exp.diameter[valid_ind]**3 / np.nansum(exp.diameter[valid_ind]**0))**(1./3),
        'D32 (um)'    : np.nansum(exp.diameter[valid_ind]**3 / np.nansum(exp.diameter[valid_ind]**2))**(1./1),
        'D43 (um)'    : np.nansum(exp.diameter[valid_ind]**4 / np.nansum(exp.diameter[valid_ind]**3))**(1./1),    
        }
    
    if use_scaled_diameters:
        if hasattr(exp,'scaled_diameter') == False:
            raise Exception('Scaled diameters array does not exist. Perform scaling first')
        
        characteristics_dict['D10 scaled (um)'] = np.nanmean(exp.scaled_diameter[valid_ind])
        characteristics_dict['D20 scaled (um)'] = np.nansum(exp.scaled_diameter[valid_ind]**2 / np.nansum(exp.scaled_diameter[valid_ind]**0))**(1./2)
        characteristics_dict['D30 scaled (um)'] = np.nansum(exp.scaled_diameter[valid_ind]**3 / np.nansum(exp.scaled_diameter[valid_ind]**0))**(1./3)
        characteristics_dict['D32 scaled (um)'] = np.nansum(exp.scaled_diameter[valid_ind]**3 / np.nansum(exp.scaled_diameter[valid_ind]**2))**(1./1)
        characteristics_dict['D43 scaled (um)'] = np.nansum(exp.scaled_diameter[valid_ind]**4 / np.nansum(exp.scaled_diameter[valid_ind]**3))**(1./1)

    return characteristics_dict

def process_full_spray_old(experiments, params, use_scaled_diameter=False):
    """
    Calculate full spray characteristics. Requires measurements made across entire cross-section of spray (r=0:spray_diameter)
    Algorithm is derived from Peter's JFM paper with modifications from Kee On. 
    Returns full spray PDF, VFDF (volume flux denstity function), and diameter characteristics

    """
    
    slot = params['slot_width']
    max_diameter = params['max_diameter']
    median_particle_diameter = params["median_particle_diameter"]
    num_dbin = params["number of diameter bins"]

    
    measurement_area = slot*1e-3 # initial guess of the measurement area

    ### Define spray radius and annulus area
    radial_pos, annulus_area, total_spray_area = compute_annulus_area(experiments)

 
    bins = np.linspace(0, max_diameter, num_dbin)
    
    num_particles_bin = np.zeros(len(bins))
    trav_len_particles_bin = np.zeros(len(bins))
    trav_time_bin = np.zeros(len(bins))
    Lm = 0
    
    ct = 0;
    dm = 0; #total diameteres
    nd = 0; #total number of particles
    ttot = 0; #total acquisition time
    nnn = 0; #number of particles at median size
    ln = 0; #traverse length at median size
    vxm = 0; #velocities at median size
    txm = 0; #gate times at median size
    
    for count, exp in enumerate(experiments):
        
        valid_int = np.where(exp.valid == 1)[0]
        
        lx = abs(exp.gate_time[valid_int] * 1e-6 * exp.velocity_chan1[valid_int]) 
    
        T = exp.time_chan1[valid_int][-1]
        if use_scaled_diameter:
            d = exp.scaled_diameter[valid_int]
        else:
            d = exp.diameter[valid_int]
            
        vx = exp.velocity_chan1[valid_int]    
        tx = exp.gate_time[valid_int] * 1e-6
    
        # finding specifically for 95 and 115 um (m - median?)    
        ind = np.where((d > median_particle_diameter[0]) & (d < median_particle_diameter[1]))
        
        if len(ind)==0:
            print("Warning, not particles found")
        else:
            ln = ln + np.nansum(vx[ind] * tx[ind])
            vxm = vxm + np.nansum(vx)
            txm = txm + np.nansum(tx)
            nnn = nnn + len(ind)
        
        tx = exp.gate_time[valid_int] * 1e-6
        
        bin_count, bin_edge = np.histogram(d, bins)
        digitized_diameter = np.digitize(d, bins)
        
        for j in range(num_dbin - 1):
            bin_ind = np.where(digitized_diameter==j+1)[0]
            num_particles_bin[j] = num_particles_bin[j] + len(bin_ind) * annulus_area[count] / measurement_area
            trav_len_particles_bin[j] = trav_len_particles_bin[j] + np.nansum(lx[bin_ind]) * annulus_area[count] / measurement_area
            trav_time_bin[j] = trav_time_bin[j] +  np.nansum(tx[bin_ind]) * annulus_area[count] / measurement_area
            
            dm = dm + np.nansum(d[bin_ind]) * annulus_area[count] / measurement_area
            Lm = Lm + np.nansum(lx[bin_ind]) * annulus_area[count] / measurement_area
            nd = nd + len(bin_ind) * annulus_area[count] / measurement_area
            #ttot = ttot + T
            
    dm = dm / nd
    LN = ln / nnn
    Lm = Lm/np.nansum(num_particles_bin) 
    trav_len_particles_bin = trav_len_particles_bin / num_particles_bin
    trav_time_bin = trav_time_bin / num_particles_bin
    
    ### Calculate F parameter - beam length adjustment
    bin_centers = bins + (bins[1] - bins[0]) / 2.
    bin_centers = bin_centers[0:-1]
    
    L1 = lambda bin_centers: (30 + 2.25 * bin_centers / (1 + 0.15 * np.exp(0.005*bin_centers)))
    
    ind = np.where(abs(bin_centers - dm) == min(abs(bin_centers - dm)))[0]
    f = (L1(bin_centers) - L1(bin_centers[ind])) / L1(bin_centers[ind])

    # compute spray parameters
    count = 0
    slit = 150e-6 
    A = Lm * slot #slit * np.sin(receiver_angle * np.pi / 180) * 750. / 250.
    V = A * Lm
    
    VFD = np.zeros((len(annulus_area), len(bin_centers)))
    VF = np.zeros((len(annulus_area), len(bin_centers)))
    NNft = np.zeros(len(bins)-1)
    
    for count, exp in enumerate(experiments):
        
        valid_int = np.where(exp.valid == 1)[0]
        if use_scaled_diameter:
            d = exp.scaled_diameter[valid_int]
        else:
            d = exp.diameter[valid_int]
        T = max(exp.time_chan1[valid_int])  # Not sure why this is end and not max
        tx = exp.gate_time[valid_int] * 1e-6
        
    
        bin_count, bin_edge = np.histogram(d, bins)
        digitized_diameter = np.digitize(d, bins)
        
        for d_bin in range(num_dbin - 1):
            bin_ind = np.where(digitized_diameter==d_bin+1)[0]
            
            print(bin_ind, d_bin)
            if len(bin_ind) == 0:
                tx_bin_mean = 0
            else:
                tx_bin_mean = trav_time_bin[d_bin]
            
            VFD[count, d_bin] = np.pi / 6 * (len(bin_ind) * (bins[d_bin]*1e-6)**3) / T / A / (1+f[d_bin]) 
            NNft[d_bin] = NNft[d_bin] + len(bin_ind) / T / A / (1+f[d_bin])*annulus_area[count] / total_spray_area # area might be wrong        
            
            VF[count, d_bin] = len(bin_ind) * (bins[d_bin]*1e-6)**3 * tx_bin_mean / T / V / (1.+f[d_bin])
    


    ### Compute Metrics to return
    
    pdfN = NNft / np.trapz(NNft, x=bin_centers)
    bin_volume = NNft * np.pi / 6 * (bins[0:-1]*1e-6)**3
    D10T=np.trapz(bin_centers**1*pdfN, x=bin_centers)**(1./1);
    D30T=np.trapz(bin_centers**3*pdfN, x=bin_centers)**(1./3);
    D32T=np.trapz(bin_centers**3*pdfN, x=bin_centers)/np.trapz(bin_centers**2*pdfN, x=bin_centers);
    D43T=np.trapz(bin_centers**4*pdfN, x=bin_centers)/np.trapz(bin_centers**3*pdfN, x=bin_centers);

    diam_dict = {"D10_full": D10T,
                 "D30_full": D30T,
                 "D32_full": D32T,
                 "D43_full": D43T};


    ### integrate VFD and VF across annulus
        
    VF_annulus = np.zeros(len(bins)-1)
    VFD_annulus = np.zeros(len(bins)-1)
    for bin_count in range(num_dbin - 1):
        VF_annulus[bin_count] = np.sum(VF[:, bin_count] * annulus_area)/total_spray_area
        VFD_annulus[bin_count] = np.sum(VFD[:, bin_count] * annulus_area)/total_spray_area
    print(dm, Lm, f)
    return pdfN, VF_annulus, VFD_annulus, VFD, bin_volume, total_spray_area, bins, diam_dict





def calc_beam_length_function(experiments, params, use_scaled_diameter=False, visualize=True):

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
    for count, exp in enumerate(experiments):
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
    
    # fit curve between beam length and diameters
    def beam_length_func(x, a, b, c):
        return a + b * x / (1 + c * np.exp(0.005 * x))
    
    popt, pcov = curve_fit(beam_length_func, bin_centers[non_zero_bins], bin_beam_length[non_zero_bins], bounds=(0, [np.inf, np.inf, np.inf]))
    beam_length_fit = beam_length_func(bin_centers, popt[0], popt[1], popt[2])
    print(popt)
    if visualize:
        fig = go.Figure()
        fig.add_trace(go.Scatter(y=bin_beam_length[non_zero_bins], x=bin_centers[non_zero_bins], mode='markers'))
        fig.add_trace(go.Scatter(y=beam_length_fit, x=bin_centers, mode='markers+lines'))
        fig.update_yaxes(title='Beam crossing length (m)')
        fig.update_xaxes(title='Diameter (um)')
        
        fig.show()
           
    # calculate the f function
    
    L1 = lambda bin_centers: (popt[0] + popt[1] * bin_centers / (1 + popt[2] * np.exp(0.005*bin_centers)))
    ind = np.where(abs(bin_centers - average_diameter) == min(abs(bin_centers - average_diameter)))[0]
    f = (L1(bin_centers) - L1(bin_centers[ind])) / L1(bin_centers[ind])
    
    return f, average_beam_length, bin_beam_length, bin_travel_time



def process_full_spray(experiments, params, use_scaled_diameter=False):
    """
    Calculate full spray characteristics. Requires measurements made across entire cross-section of spray (r=0:spray_diameter)
    Algorithm is derived from Peter's JFM paper with modifications from Kee On. 
    Returns full spray PDF, VFDF (volume flux denstity function), and diameter characteristics
    
    I've removed the divison by the full spray area. We should just be integrating the VFD across the spray. that gives us volume/s for the full spray. 
    No need to divide again by the total spray area.

    """
    
    measurement_width = params['slot_width']
    max_diameter = params['max_diameter']
    median_particle_diameter = params["median_particle_diameter"]
    num_dbin = params["number of diameter bins"]
    average_beam_length = params['average beam length']
    bin_travel_time = params['bin travel time']
    f = params['beam length function']

    bins = np.linspace(0, max_diameter, num_dbin)
    bin_centers = bins + (bins[1] - bins[0]) / 2.
    bin_centers = bin_centers[0:-1]

    ### Define spray radius and annulus area
    radial_pos, annulus_area, total_spray_area = compute_annulus_area(experiments)

 
    # compute spray parameters
    count = 0
    A = average_beam_length * measurement_width #slit * np.sin(receiver_angle * np.pi / 180) * 750. / 250.
    V = A * average_beam_length
    
    VFD = np.zeros((len(annulus_area), len(bin_centers)))
    VF = np.zeros((len(annulus_area), len(bin_centers)))
    NFD = np.zeros(len(bins)-1)
    
    for count, exp in enumerate(experiments):
        
        valid_int = np.where(exp.valid == 1)[0]
        if use_scaled_diameter:
            d = exp.scaled_diameter[valid_int]
        else:
            d = exp.diameter[valid_int]
        T = max(exp.time_chan1[valid_int])  # Not sure why this is end and not max
        tx = exp.gate_time[valid_int] * 1e-6
        
    
        bin_count, bin_edge = np.histogram(d, bins)
        digitized_diameter = np.digitize(d, bins)
        
        for d_bin in range(num_dbin - 1):
            bin_ind = np.where(digitized_diameter==d_bin+1)[0]
            
            print(bin_ind, d_bin)
            if len(bin_ind) == 0:
                tx_bin_mean = 0
                bin_particle_count = 0
            else:
                tx_bin_mean = bin_travel_time[d_bin]
                bin_particle_count = len(bin_ind)
            
            VFD[count, d_bin] = np.pi / 6 * (bin_particle_count * (bins[d_bin]*1e-6)**3) / T / A / (1+f[d_bin]) 
            NFD[d_bin] = NFD[d_bin] + bin_particle_count / T / A / (1+f[d_bin])*annulus_area[count] # area might be wrong        
            
            VF[count, d_bin] = bin_particle_count * (bins[d_bin]*1e-6)**3 * tx_bin_mean / T / V / (1.+f[d_bin])
    



    ### Compute Metrics to return
    
    pdfN = NFD / np.trapz(NFD, x=bin_centers)
    bin_volume = NFD * np.pi / 6 * (bins[0:-1]*1e-6)**3
    D10T=np.trapz(bin_centers**1*pdfN, x=bin_centers)**(1./1);
    D30T=np.trapz(bin_centers**3*pdfN, x=bin_centers)**(1./3);
    D32T=np.trapz(bin_centers**3*pdfN, x=bin_centers)/np.trapz(bin_centers**2*pdfN, x=bin_centers);
    D43T=np.trapz(bin_centers**4*pdfN, x=bin_centers)/np.trapz(bin_centers**3*pdfN, x=bin_centers);

    diam_dict = {"D10_full": D10T,
                 "D30_full": D30T,
                 "D32_full": D32T,
                 "D43_full": D43T};


    ### integrate VFD and VF across annulus
        
    VF_annulus = np.zeros(len(bins)-1)
    VFD_annulus = np.zeros(len(bins)-1)
    for bin_count in range(num_dbin - 1):
        VF_annulus[bin_count] = np.sum(VF[:, bin_count] * annulus_area)
        VFD_annulus[bin_count] = np.sum(VFD[:, bin_count] * annulus_area)
    print(T)
    return pdfN, VF_annulus, VFD_annulus, VFD, bin_volume, total_spray_area, bins, diam_dict


