#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Mon 15 Feb 2021
@author: goharshoukat
'''

import math
import numpy as np
import Linear_Airy_Wave_Solution
import scipy.interpolate as interpolate

class Morisons(Linear_Airy_Wave_Solution.Set_Wave_Field):
    
    def __init__(self, depth, time_period, tide_velocity, amplitude, diameter):
        super().__init__(depth, time_period, tide_velocity, amplitude)
        self.diameter = diameter
        
    
    def elevation(self, x: float, time: np.ndarray):
        if type(x) == np.ndarray and type(time) == np.ndarray:
            raise TypeError ('x can be a float value only')
        else:
            return   (self.amplitude * np.cos(self.kfromw()\
                                                * x - self.omega*time))  
    
    #Determines if Morison's Equation is applicable or not. 
    #Makes use of the displacement parameter for judgement. 
    #KC can also be used. KC = 2pi A/D
    def diffraction_effects(self):
        ratio = self.diameter/self.wave_length()
        if ratio < 0.2:
            print('Small Body Assumption')
        else:
            print('Large Body Assumption')
        return None
    
    
    def z_array(self, total_depth, delta_z, start =0):
        return np.arange(start+delta_z, -total_depth, -delta_z).reshape(-1,1)
    
    def time(self, total_time, delta_t, start = 0):
        return np.arange(start, total_time+delta_t, delta_t).reshape(-1,1)
    
    #time is passed on as a numpy array of n x 1 dimesnions
    #z is passed on a numpy arraym x 1 dimensions
    #z goes from 0 to depth. 
    def horizontal_velocity(self, x, time, z):
        self.U = self.amplitude * self.omega * np.dot((np.cosh(self.kfromw()\
                * (z + self.depth)) / np.sinh(self.kfromw() * self.depth)), \
                            (np.cos(self.kfromw() * x - self.omega * time)).T)
        return self.U
    #different function from the one in the linear airy wave solutions library because in that
    #the functions are spatially discretised. Here they are discretised in time. 
    def horizontal_acceleration(self, x, time, z):
        self.accel = self.amplitude * pow(self.omega,2) *\
            np.dot((np.cosh(self.kfromw() * (z + self.depth))/\
            np.sinh(self.kfromw() * self.depth)), \
                   (np.sin(self.kfromw() * x - self.omega * time)).T)
        return self.accel
    
    def Reynolds(self):
        return (self.U * self.diameter / self.nu)
    
    def kc(self):
        return (self.U * self.time_period / self.diameter) 
    
    #array can not be passed on to this function. loop it through
    def coefficient(self, local_Re, local_kc):
       # shape_local_Re = np.shape(local_Re)
       # shape_local_kc = np.shape(local_kc)
     #   local_Re = local_Re.flatten('F')
      #  local_kc = local_kc.flatten('F')
        if type(local_Re) == np.ndarray or type(local_kc) == np.ndarray:
            raise TypeError ('x can be a non-array value only') 
        Re = np.array((-10, 1e1, 1e3, 5e5, 1e10))
        kc = np.array((-0.05, 5, 15, 50, 100))
        KC, RE = np.meshgrid(kc, Re)
        Cm = np.array([[2, 2, 0.6, 1, 1.1], 
                       [2, 2, 0.6, 1, 1.2],
                       [2, 2, 0.8, 1.4, 1.5],
                       [2, 2, 1.3, 1.6, 1.7], 
                       [2, 2, 2, 2, 2]])
        
        Cd = np.array([[2, 2, 2.5, 1.5, 1.3],
                       [2, 2, 2.5, 1.5, 1.3],
                       [2, 2, 2.5, 1.5, 1.3],
                       [2, 1.3, 1, 0.6, 0.5],
                       [2, 0.5, 0.7, 0.6, 0.6]])
        
        #value isnt exactly accurate after interpolation. this function has some errors from scipy
        ip_Cm = interpolate.interp2d(KC, RE, Cm, kind = 'cubic')
        local_Cm = ip_Cm(local_kc, local_Re)
        
        ip_Cd = interpolate.interp2d(KC, RE, Cd, kind = 'cubic')
        local_Cd = ip_Cd(local_kc, local_Re)
       # local_Re = local_Re.reshape(shape_local_Re[0], shape_local_Re[1])
      #  local_kc = local_kc.reshape(shape_local_kc[0], shape_local_kc[1])
        return local_Cm, local_Cd
    
    def force(self, local_Cd, local_Cm):
        f_drag = 0.5 * self.rho * local_Cd.T * self.diameter * abs(self.U) * self.U
        f_inertial = (self.rho * local_Cm.T * math.pi * pow(self.diameter),2)\
            * self.accel
            
        f_total = f_drag + f_inertial
        force = {f_drag: 'Drag Force', f_inertial: 'Inertial Force', f_total: 'Total Force'}
        return force
        
    
field = Morisons(50, 10, 0, 0.3, 10)
t = field.time(20,0.5)
x = 0
elevation = field.elevation(0,t)
field.diffraction_effects()
z_array = field.z_array(50, 0.5)
hor_velocity = field.horizontal_velocity(x, t, z_array)
a = field.horizontal_acceleration(x, t, z_array)
Re = field.Reynolds()
kc = field.kc()
Cm, Cd = field.coefficient(10000, 10)


