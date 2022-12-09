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
    
    
    def __init__(self, name, params, data_dir):
        
        self._name = name
        self._path = data_dir + params['Data folder name']
        self.tidal_volume = params['Tidal volume (mL)']
        self.inspiratory_period = params["Inspiration period (s)"]
        self.measurement_distance = params['Centerline distance (mm)']
        self.droplet_diameter_scale = 0.467
        self.minintensity = params['minintensity']
        self.maxdiamdiff = params['maxdiamdiff']
        self.a1 = params['a1']
        self.c1 = params['c1']
        self.ratio = params['ratio']
        self.c2 = params['c2']
 
        self.estimate_flowrate()
        
        self.load_pdpa_data()
        
    def load_pdpa_data(self):
        
        data = pd.read_csv(self._path, delimiter=',', skiprows=1)
        
        # self.time_chan1 = data['Time Ch. 1 (sec)']
        # self.time_chan2 = data[' "Time Ch. 2 (sec)"']
        self.velocity_chan1 = np.array(data[' "Velocity Ch. 1 (m/sec)"'].values * -1)
        # self.velocity_chan2 = np.array(data[' "Velocity Ch. 2 (m/sec)"'])
        self.intensity = np.array(data[' "Intensity (mV)"'])
        self.diameter_diff = np.array(data[' "Diameter Dif (um)"'])
        
        diameter = data[' "Diameter (um)"'].values
        self.diameter = np.array([float(x) if x!=' Fillin' else np.NaN for x in diameter])
        
        self.volume = 4. / 3. * np.pi * (self.diameter / 2.)**3
        
        zero_intensity_count = len(np.where(self.intensity == 0)[0])
        print("Percentage of zero intensity samples: {:3f}%".format( zero_intensity_count / len(self.intensity) * 100 ))
        
    def scale_diameter(self):
        print("Scaling diameters by a factor of:", self.droplet_diameter_scale)
        self.diameter = self.diameter * self.droplet_diameter_scale
    
    def estimate_flowrate(self):
        self.estimated_flowrate = self.tidal_volume / self.inspiratory_period * 60 / 1000
        
        
    def filter_diameters(self, filter_curves=True):
        
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
            fig.update_layout(title="{:s}; a1={:0.2f}; c1={:0.2f}; ratio={:0.2f}; c2={:0.2f}".format(self._name, self.a1, self.c1, self.ratio, self.c2))
            
            # add a save image
            
            fig.show()
            
        valid_final = (valid1 & valid0) == True
        idx_valid = np.nonzero(valid_final)[0]
        
        percent_valid = 100*len(idx_valid)/len(self.diameter)
        print("{:2.2f} % valid particles".format(percent_valid))
        
        
        self.diameter = self.diameter[valid_final]
        self.velocity_chan1 = self.velocity_chan1[valid_final]
        self.intensity = self.intensity[valid_final]
        
        
        print(len(self.diameter))
        
    def remove_zero_intensity(self):
        
        zero_intensity_values = np.where(self.intensity == 0)
        self.diameter = np.delete(self.diameter, zero_intensity_values)
        self.velocity_chan1 = np.delete(self.velocity_chan1, zero_intensity_values)
        self.intensity = np.delete(self.intensity, zero_intensity_values)
        
        
        
        

        
        
        
        
        
        