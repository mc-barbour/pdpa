#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 09:52:14 2022

@author: mbarbour
"""

import pandas as pd
import numpy as np

import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots


pio.renderers.default = "browser"


class Run():
    
    
    def __init__(self, name, params, data_dir, raw_bool=True):
        
        self._name = name
        self._path = data_dir + params['Data folder name']
        self.data_dir = data_dir
        self.raw_bool = raw_bool
        
        self.tidal_volume = params['Tidal volume (mL)']
        self.inspiratory_period = params["Inspiration period (s)"]
        self.measurement_distance = params['Centerline distance (mm)']
        self.droplet_diameter_scale = 0.467
        self.minintensity = params['minintensity']
        self.maxdiamdiff = params['maxdiamdiff']
        self.maxdiameter = params['maxdiameter']
        self.a1 = params['a1']
        self.c1 = params['c1']
        self.ratio = params['ratio']
        self.c2 = params['c2']
 
        self.estimate_flowrate()
        
        self.load_pdpa_data()
        
        self.centerline_distance = params['Centerline distance (mm)']
        self.radial_position = params['Radial position (mm)']
        
    def load_pdpa_data(self):
        
        if self.raw_bool:
            data = pd.read_csv(self._path, delimiter=',', skiprows=1)
        else:
            data = pd.read_csv(self._path, delimiter=',', skiprows=0)
            self.valid = data[' "Validity"'].values
        
        if ' "Scaled Diameter"' in data.columns:
            print('Loading scaled diameters')
            self.scaled_diameter = data[' "Scaled Diameter"'].values
        
        self.time_chan1 = data['Time Ch. 1 (sec)'].values
        self.time_chan2 = data[' "Time Ch. 2 (sec)"'].values
        self.gate_time = data[' "Gate Time Ch. 1 (usec)"'].values
        
        self.velocity_chan1 = np.array(data[' "Velocity Ch. 1 (m/sec)"'].values)
        self.velocity_chan2 = np.array(data[' "Velocity Ch. 2 (m/sec)"'])
        self.intensity = np.array(data[' "Intensity (mV)"'])
        self.diameter_diff = np.array(data[' "Diameter Dif (um)"'])
        
        diameter = data[' "Diameter (um)"'].values
        self.diameter = np.array([float(x) if x!=' Fillin' else np.NaN for x in diameter])
        
        self.volume = 4. / 3. * np.pi * (self.diameter / 2.)**3
        
        zero_intensity_count = len(np.where(self.intensity == 0)[0])
        print("Percentage of zero intensity samples: {:3f}%".format( zero_intensity_count / len(self.intensity) * 100 ))
        
    def scale_diameter(self):
        print("Scaling diameters by a factor of:", self.droplet_diameter_scale)
        self.scaled_diameter = self.diameter * self.droplet_diameter_scale
    
    def estimate_flowrate(self):
        self.estimated_flowrate = self.tidal_volume / self.inspiratory_period * 60 / 1000
        
        
    def filter_diameters(self, filter_curves=True, max_diameter=True, show_plot=True):
        
        # remove nans
        valid0 = ~np.isnan(self.diameter)
        
        if filter_curves:
            print(max(self.diameter))
            x = np.arange(0, np.ceil(np.nanmax(self.diameter)), 1)
            
            y1 = self.a1 * x**2 + self.c1
            i1 = self.a1 * self.diameter**2 + self.c1
            
            y2 = self.ratio * self.a1 * x**2 + self.c2
            i2 = self.ratio * self.a1 * self.diameter**2 + self.c2
        
            valid1 = (self.intensity < i1) & (self.intensity > i2) == True
            print(len(i1), len(i2), len(self.intensity))
            not_valid = valid1==False
            
            if show_plot:
                fig = go.Figure()
                fig.add_trace(go.Scatter(x=self.diameter[valid1], y=self.intensity[valid1],
                                         mode='markers', marker_color='blue',
                                         name='valid'))
                fig.add_trace(go.Scatter(x=self.diameter[not_valid], y=self.intensity[not_valid],
                                         mode='markers', marker_color='red',
                                         name='not valid'))
                fig.add_trace(go.Scatter(x=x, y=y1, name='Curve 1', marker_color='black', line_dash='dash'))
                fig.add_trace(go.Scatter(x=x, y=y2, name='Curve 2',marker_color='black', line_dash='dash'))
                fig.update_yaxes(title='Intensity', range=[self.c2,max(self.intensity)])
                fig.update_layout(title="{:s}; a1={:0.3f}; c1={:0.3f}; ratio={:0.3f}; c2={:0.2f}".format(self._name, self.a1, self.c1, self.ratio, self.c2))
                
                # add a save image
                
                fig.show()
        if max_diameter:
            valid2 = (self.diameter < self.maxdiameter)
            valid_final = (valid1 & valid0 & valid2) == True
        else:
            valid_final = (valid1 & valid0) == True
            
        idx_valid = np.nonzero(valid_final)[0]
        
        percent_valid = 100*len(idx_valid)/len(self.diameter)
        print("{:2.2f} % valid particles".format(percent_valid))
        
        self.valid = valid_final * 1 #convert from boolean to [0,1]
        print(valid_final)
        # self.diameter = self.diameter[valid_final]
        # self.velocity_chan1 = self.velocity_chan1[valid_final]
        # self.intensity = self.intensity[valid_final]
        
        
        print(len(self.diameter))
        
    def remove_zero_intensity(self):
        
        zero_intensity_values = np.where(self.intensity == 0)
        self.diameter = np.delete(self.diameter, zero_intensity_values)
        self.velocity_chan1 = np.delete(self.velocity_chan1, zero_intensity_values)
        self.intensity = np.delete(self.intensity, zero_intensity_values)
        
    def filter_maxdiameter(self):
        valid = (self.diameter < self.maxdiameter)
        idx_valid = np.nonzero(valid)[0]
        # if valid in dir(self):
        #     self.valid.append(valid)
        
        percent_valid = 100*len(idx_valid)/len(self.diameter)
        print("{:2.2f} % valid particles".format(percent_valid))
        self.valid = valid*1
        
    def save_filtered_run(self):
       
        data = pd.read_csv(self._path, delimiter=',', skiprows=1)
        
        diameter = data[' "Diameter (um)"'].values
        self.diameter = np.array([float(x) if x!=' Fillin' else np.NaN for x in diameter])
        
        data[' "Diameter (um)"'] = self.diameter
        data[' "Validity"'] = self.valid
        
        if hasattr(self, 'scaled_diameter'):
            data[' "Scaled Diameter"'] = self.scaled_diameter
        else:
            print("No scaled diameter to save")
            
        process_data_dir = self.data_dir + "../valid_data/"
        filename = self._path.split("/")[-1]
        
        print("Saving valid particles to: ", process_data_dir + filename)
        
        data.to_csv(process_data_dir + filename)
        
        
        
        
        
        
        

        
        
        
        
        
        