# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:44:12 2013

@author: eric

This python script is imported to microcavity_TM_approach_N_layers.py, and provides the code which
imports the data from the config file

Ensure that the "configobj" package is installed on your computer
"""
from configobj import ConfigObj
from validate import Validator
import sys


config_file = 'config.ini'
configspec_file = 'configspec.ini'


# read in configuration file and validate against configspec
config = ConfigObj(config_file, configspec=configspec_file)
validator = Validator()
if config.validate(validator) != True:
    print "Error reading config file - input parameters not valid."
    sys.exit(1)
#-------------------------------------------------------------#
#read run parameters
polarization = config['run_params']['polarization']
inc_angle = config['run_params']['inc_angle'] #in degrees
lambda_0 = config['run_params']['lambda_0'] 
lambda_start = config['run_params']['lambda_start'] 
lambda_end = config['run_params']['lambda_end'] 
num_points = config['run_params']['num_points']

#-------------------------------------------------------------#
#read sample specification
materials = config['sample']['available_materials']
refractive_indices = config['sample']['available_indices_refraction']
substrate = config['sample']['substrate']
medium = config['sample']['medium']
layer_list = config['sample']['layers']
thickness_list = config['sample']['thickness']

#-------------------------------------------------------------#
#Input Validation

if len(layer_list) != len(thickness_list):
    print "Error - layers list and thickness list do not have the same length."
    sys.exit(1)

if len(materials) != len(refractive_indices):
    print "Error - available materials list and available indices of refraction list do not have the same length."
    sys.exit(1)

for i in range(len(thickness_list)):
    if thickness_list[i]<=0:
        print "Error - thicknesses must be non-negative values"
        sys.exit(1)

if lambda_end <= lambda_start:
    print "Error - the end-wavelength must be greater than the start-wavelength"
    sys.exit(1)

if lambda_0 < lambda_start or lambda_0 > lambda_end:
    print "Warning - the design wavelength falls outside of the range of the scan-wavelengths\nContinuing..."

#-------------------------------------------------------------#
# Add parameters into dicitonaries
material_dictionary = dict(zip(materials, refractive_indices))
param_dictionary = {'polarization':polarization, 'inc_angle':inc_angle, 'lambda_0':lambda_0, 'lambda_start':lambda_start, 'lambda_end':lambda_end, 'num_points':num_points, 'substrate':substrate, 'medium':medium, 'layer_list':layer_list, 'thickness_list':thickness_list}