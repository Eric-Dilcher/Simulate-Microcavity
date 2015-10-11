# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:44:12 2013

@author: eric

This python script is imported to microcavity_TM_approach_N_layers.py, and provides the code which
imports the data from the config file

Ensure that the "configobj" package is installed on your computer
"""
from math import pi
from configobj import ConfigObj
from validate import Validator
import sys


config_file = 'config.ini'#enter the filename of the config file
configspec_file = 'configspec.ini'#enter the filename of the configspec file


# read in configuration file and validate against configspec
config = ConfigObj(config_file, configspec=configspec_file)
validator = Validator()
if config.validate(validator) != True:
    print "Error reading config file"
    sys.exit(1)
#-------------------------------------------------------------#
#read run_params
polarization = config['run_params']['polarization']
inc_angle = config['run_params']['inc_angle'] #incident angle in degrees
lambda_0 = config['run_params']['lambda_0'] #design wavelength of Bragg stack
lambda_start = config['run_params']['lambda_start'] #start wavelength
lambda_end = config['run_params']['lambda_end'] #end wavelength
num_points = config['run_params']['num_points'] #number of samples

#-------------------------------------------------------------#
#read sample specification
materials = config['sample']['available_materials']
refractive_indices = config['sample']['available_indices_refraction']
substrate = config['sample']['substrate']#string that holds the name of the substrate
medium = config['sample']['medium']#string that holds the name of the medium surrounding the sample
layer_list = config['sample']['layers'] #string_list that holds the names of the layers
thickness_list = config['sample']['thickness'] #float_list that holds the thicknesses of the layers (in fractions of design wavelengths)

#Input Validation
#check to make sure layer name list, and layer thickness list have te same length
if len(layer_list) != len(thickness_list):
    print "Error - layers list and thickness list do not have the same number of items"
    sys.exit(1)
    
#Input Validation
#check to make sure layer name list, and layer thickness list have te same length
if len(materials) != len(refractive_indices):
    print "Error - available materials list and available indices of refraction do not have the same number of items."
    sys.exit(1)
#check to make sure that the thickness values are non-negative
for i in range(len(thickness_list)):
    if thickness_list[i]<=0:
        print "Error - thicknesses must be non-negative values"
        sys.exit(1)
        
#check to make sure that the end-wavelength is greater than the start wavelength
if lambda_end < lambda_start:
    print "Error - the end-wavelength must be greater than the start-wavelength"
    sys.exit(1)
    
#give a warning if the design wavelength falls outside of the range defined by the start and end wavelengths
if lambda_0 < lambda_start or lambda_0 > lambda_end:
    print "Warning - the design wavelength falls outside of the range of the scan-wavelengths\nContinuing"

material_dictionary = dict(zip(materials, refractive_indices))
param_dictionary = {'polarization':polarization, 'inc_angle':inc_angle, 'lambda_0':lambda_0, 'lambda_start':lambda_start, 'lambda_end':lambda_end, 'num_points':num_points, 'substrate':substrate, 'medium':medium, 'layer_list':layer_list, 'thickness_list':thickness_list}