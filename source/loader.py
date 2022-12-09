#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 10:09:04 2022

@author: mbarbour
"""
import os
import pandas as pd

from source.experiment import Run

def get_experiment_names(data_dir):
    
    """"Print all available experiments in the experimental log"""
    
    experiment_log_filename = data_dir + 'experiment_log.xlsx'
    
    if os.path.exists(experiment_log_filename):
        experiment_log = pd.read_excel(data_dir + 'experiment_log.xlsx', index_col='Name')

    else:
        raise Exception("Error: No experimental log found in director: ", data_dir)
        
    names = []
    for name in experiment_log.index.values:
        names.append(name)
        print("{:s}".format(name))
    return names


def load_experiment(name, data_dir):
    
    experiment_log_filename = data_dir + 'experiment_log.xlsx'
    
    if os.path.exists(experiment_log_filename):
        experiment_log = pd.read_excel(data_dir + 'experiment_log.xlsx', index_col='Name')
    
    else:
        raise Exception("Error: No experimental log found in director: ", data_dir)
        
        
    if name in experiment_log.index.values:
        print(experiment_log.loc[name])
        return Run(name, experiment_log.loc[name], data_dir)
    else:
        raise ValueError("No experiment found with name {:s}. Check spelling and paths in experiment.py".format(name))