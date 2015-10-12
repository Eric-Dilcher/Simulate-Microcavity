from scipy.constants import c, pi
from scipy import matrix, sin, cos, arcsin
import read_config as init

RAD_PER_DEG = pi/180

material_dictionary = init.material_dictionary
#class to define a layer in a multilayer stack
class Layer(object):
    def __init__(self, material, thick_frac, param_dict):#Layer class constructor
        #material name
        self.material_name = material
        #thickness of layer as a fraction of design wavelength
        self.thick_frac = thick_frac
        #index of refraction
        self.index_r = material_dictionary[material]
        self.set_thick(param_dict)
        self.set_theta(param_dict)
        self.set_gamma_perp()
        self.set_gamma_para()
        
    #method to calculate and the layer thickness
    def set_thick(self, param_dict):
        self.thick = self.thick_frac*param_dict['lambda_0']/self.index_r
        
    #method to calculate and set the angle of outgoing light using Snell's law
    def set_theta(self, param_dict):
        self.theta = arcsin((material_dictionary[param_dict['medium']]/self.index_r)*sin(RAD_PER_DEG*param_dict['inc_angle']))
        
    #method to calculateand set delta
    def set_delta(self, wlength):
        if self.thick_frac <= 0:
            self.delta = None
        else:
            self.delta = 2*pi*self.index_r*self.thick*cos(self.theta)/wlength
        
    #method to find gamma factors
    def set_gamma_perp(self):
        self.gamma_perp = self.index_r*cos(self.theta)/c
        
    def set_gamma_para(self):
        self.gamma_para = self.index_r/(cos(self.theta)*c)
    #method to set the transfer matrices
    def set_tmatrix_perp(self):
        if self.delta == None:
            self.tmatrix_perp = None
        else:
            self.tmatrix_perp = matrix([[cos(self.delta),1j*sin(self.delta)/self.gamma_perp], 
                                    [1j*self.gamma_perp*sin(self.delta), cos(self.delta)]])
    def set_tmatrix_para(self):
        if self.delta == None:
            self.tmatrix_para = None
        else:
            self.tmatrix_para = matrix([[cos(self.delta),1j*sin(self.delta)/self.gamma_para], 
                                    [1j*self.gamma_para*sin(self.delta), cos(self.delta)]])

    #method to refresh delta values, and transfer matrix of layer object    
    def refresh_values(self, wlength, param_dict):
        #redefine delta taking into account a different incident wavelength
        self.set_delta(wlength)
        #refresh the appropriate tmatrices
        if param_dict['polarization']==-1:#if the light is unpolarized
            self.set_tmatrix_perp()
            self.set_tmatrix_para()
        elif param_dict['polarization']==0:#if the light is polarized perpendicular to the plane of incidence
            self.set_tmatrix_perp()
        elif param_dict['polarization']==1:#if the light is polarized parallel to the plane of incidence
            self.set_tmatrix_para()
        else:
            print "Error - polarization indicator not in range"