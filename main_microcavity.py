# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:15:06 2013

@author: Eric Dilcher

This code relies on the transfer matrix theory (see "Introduction to Optics" by Pedrotti)
to characterize the reflective and transmissive properties of multilayer films.

Any types, number of, thicknesses of films can be analyzed by entering the characteristics of the film in the config file.

There is a config file 'config.ini' where the user is to enter the parameters of the simulation.
If the user would like to expand the allowable materials, update the Material_Dictionary below.  In addition
update the comments in config.ini, and configspec.ini

Note:  the package "configobj" needs to be installed before this code can be written
"""

from scipy import matrix
import numpy as np
from matplotlib import pyplot as plt
import read_config as init
from Layer import Layer

#--------------------------------------------------------------#    
# Functions

def create_layer_objs(param_dict):
    layers_objs = []
    for layer, thick in zip(param_dict["layer_list"], param_dict["thickness_list"]):
        layers_objs.append(Layer(layer, thick, param_dict))
    return layers_objs

def calculate(param_dictionary, medium, substrate, layer_objs):
    start = param_dictionary["lambda_start"]
    end = param_dictionary["lambda_end"]
    num_points = param_dictionary["num_points"]
    polarization = param_dictionary["polarization"]
    
    lambda_values = np.linspace(start, end, num_points)
    transm_coeff = []
    refl_coeff = []

    for wavelen in lambda_values:
        update_all_layers(wavelen, param_dictionary, layer_objs)#recalculate the values taking into account the changed wavelength
        if polarization == -1:
            T_perp, R_perp = calc_trans_refl_at_wavelen_perp(wavelen, medium, substrate, layer_objs)
            T_para, R_para = calc_trans_refl_at_wavelen_para(wavelen, medium, substrate, layer_objs)
            transm_coeff.append(np.average([T_perp, T_para]))
            refl_coeff.append(np.average([R_perp, R_para]))
        elif polarization == 0:
            T, R = calc_trans_refl_at_wavelen_perp(wavelen, medium, substrate, layer_objs)
            transm_coeff.append(T)
            refl_coeff.append(R)
        else:
            T, R = calc_trans_refl_at_wavelen_para(wavelen, medium, substrate, layer_objs)
            transm_coeff.append(T)
            refl_coeff.append(R)
    
    #transform to numpy arrays
    transm_coeff = np.array(transm_coeff)
    refl_coeff = np.array(refl_coeff)
    return (lambda_values, transm_coeff, refl_coeff)

#update all the layer values taking into account the changed wavelength
def update_all_layers(wavelen, param_dict, layer_objs):
    for layer in layer_objs:
        layer.refresh_values(wavelen, param_dict)

#--------------------------------------------------------------#    
# Perpendicular polarization specific functions

def calc_trans_refl_at_wavelen_perp(wavelen, medium, substrate, layer_objs):
    sysmat = calc_sysmat_perp(layer_objs)#recalculate the system matrix
    t = calc_transmission_perp(medium, substrate, sysmat)#recalculate the transmission coefficient
    T=(np.abs(t)**2)*100#recalculate the transmission amplitude in percent
    r = calc_reflection_perp(medium, substrate, sysmat)#recalculate the reflection coefficient
    R = (np.abs(r)**2)*100#recalculate the reflection amplitude in percent
    return (T,R)
    
#calculate the system matrix
def calc_sysmat_perp(layer_objs):
    sysmat = matrix([[1.0, 0], [0, 1.0]])
    for layer in layer_objs:
        sysmat = sysmat*layer.tmatrix_perp
    return sysmat

#function to find the transmission coefficient of a stack of layers
def calc_transmission_perp(medium, substrate, sysmat_perp):
    g_0 = medium.gamma_perp
    g_s = substrate.gamma_perp
    return 2*g_0/(g_0*sysmat_perp[0,0]+g_0*g_s*sysmat_perp[0,1]+sysmat_perp[1,0]+g_s*sysmat_perp[1,1])

#function to find the reflection coefficient of a stack of layers
def calc_reflection_perp(medium, substrate, sysmat_perp):
    g_0 = medium.gamma_perp
    g_s = substrate.gamma_perp
    return (g_0*sysmat_perp[0,0]+g_0*g_s*sysmat_perp[0,1]-sysmat_perp[1,0]-g_s*sysmat_perp[1,1])/(g_0*sysmat_perp[0,0]+g_0*g_s*sysmat_perp[0,1]+sysmat_perp[1,0]+g_s*sysmat_perp[1,1])

#--------------------------------------------------------------#    
# Parallel polarization specific functions

def calc_trans_refl_at_wavelen_para(wavelen, medium, substrate, layer_objs):
    sysmat = calc_sysmat_para(layer_objs)#recalculate the system matrix
    t = calc_transmission_para(medium, substrate, sysmat)#recalculate the transmission coefficient
    T=(np.abs(t)**2)*100#recalculate the transmission amplitude in percent
    r = calc_reflection_para(medium, substrate, sysmat)#recalculate the reflection coefficient
    R = (np.abs(r)**2)*100#recalculate the reflection amplitude in percent
    return (T,R)
    
def calc_sysmat_para(layer_objs):
    sysmat = matrix([[1.0, 0], [0, 1.0]])
    for layer in layer_objs:
        sysmat = sysmat*layer.tmatrix_para
    return sysmat
    
def calc_transmission_para(medium, substrate, sysmat_para):
    g_0 = medium.gamma_para
    g_s = substrate.gamma_para
    return 2*g_0/(g_0*sysmat_para[0,0]+g_0*g_s*sysmat_para[0,1]+sysmat_para[1,0]+g_s*sysmat_para[1,1])

def calc_reflection_para(medium, substrate, sysmat_para):
    g_0 = medium.gamma_para
    g_s = substrate.gamma_para
    return (g_0*sysmat_para[0,0]+g_0*g_s*sysmat_para[0,1]-sysmat_para[1,0]-g_s*sysmat_para[1,1])/(g_0*sysmat_para[0,0]+g_0*g_s*sysmat_para[0,1]+sysmat_para[1,0]+g_s*sysmat_para[1,1])


def write_file(param_dictionary, wavelens, transmissions, reflections, filename = ""):
    polarization_dict = {-1: "Unpolarized", 0: "Perpendicular", 1: "Parallel"}
    wavelens, transmissions, reflections = np.array(wavelens), np.array(transmissions), np.array(reflections)
    data = np.c_[wavelens, transmissions, reflections] #stack arrays horizontally (3 colums)
    
    start = param_dictionary["lambda_start"]
    end = param_dictionary["lambda_end"]
    increment = (1.0*end - start)/param_dictionary["num_points"]
    center_wavelen = param_dictionary["lambda_0"]
    polarization = param_dictionary["polarization"]
    
    if filename == "":
        filename = "microcavity_output_"+str(center_wavelen)+"nm.txt"
    
    header = ("Start wavelength: "+str(start)+"\tEnd wavelength: "+str(end)+"\tIncrement by: "+str(increment)+ "\n"
        + "Design wavelength: "+str(center_wavelen) + "\tPolarization: "+ polarization_dict[polarization] + "\n"
        + "Wavelength (nm) \tTransmission (%)\tReflection (%)")
    np.savetxt(filename, data, delimiter = "\t", header = header )

def plot_data(wavelens, transm_coeff, refl_coeff, param_dictionary):
    polarization_dict = {-1: "Unpolarized", 0: "Perpendicular", 1: "Parallel"}
    center_wavelen = param_dictionary["lambda_0"]
    incident_angle = param_dictionary["inc_angle"]
    polarization = param_dictionary["polarization"]
    title = ("Design Wavelength: %d  Incident Angle: %d degrees\n" %(center_wavelen, incident_angle)
        +"Polarization: " + polarization_dict[polarization])
    plt.plot(wavelens, transm_coeff, label = "Transmission")
    plt.plot(wavelens, refl_coeff, label = "Reflection")
    plt.ylim(-5,105)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Transmission/Reflection (%)")
    plt.title(title)
    plt.legend()
    plt.show()
#-------------------------------------------------------------#
# Simulation Setup
param_dictionary = init.param_dictionary
medium = Layer(param_dictionary["medium"], -1, param_dictionary)#create a layer object for the medium 
substrate = Layer(param_dictionary["substrate"], -1, param_dictionary)#create a layer object for the substrate
layer_objs = create_layer_objs(param_dictionary)#create a list of layer objects using the create_layer_objects function

#--------------------------------------------------------------#
# Run Simulation

lambda_values, transm_coeff, refl_coeff = calculate(param_dictionary, medium, substrate, layer_objs)
write_file(param_dictionary, lambda_values, transm_coeff, refl_coeff)
plot_data(lambda_values, transm_coeff, refl_coeff, param_dictionary)


