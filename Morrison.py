#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Mon 15 Feb 2021
@author: goharshoukat
'''

import math
import numpy as np
import Linear_Airy_Wave_Solution

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
        return self.amplitude * self.omega * np.dot((np.cosh(self.kfromw() * (z + self.depth))/\
            np.sinh(self.kfromw() * self.depth)), (np.sin(self.kfromw() * x - self.omega * time)).T)
    

calc = Morisons(50, 5, 0, 0.1, 5)
z = calc.z_array(50, 0.1)
t = calc.time(50, 0.1)

