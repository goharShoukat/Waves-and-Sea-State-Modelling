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
        self.U = self.amplitude * self.omega * np.dot((np.cosh(self.kfromw() * (z + self.depth))/\
            np.sinh(self.kfromw() * self.depth)), (np.sin(self.kfromw() * x - self.omega * time)).T)
        return self.U
    
    def Reynolds(self):
        return (self.U * self.diameter / self.nu)
    
    def kc(self):
        return (2 * np.pi * self.amplitude / self.diameter) 
    
    def coefficient(self, local_Re, local_kc):
        Re = np.array((-10, 1e1, 1e3, 5e5, 1e10))
        kc = np.array(-0.05, 5, 15, 50, 100)
        KC, RE = np.meshgrid(kc, Re)
        Cm = np.array([[2, 2, 0.6, 1, 1.1], 
               [2, 2, 0.6, 1, 1.2],
               [2, 2, 0.8, 1.4, 1.5],
               [2, 2, 1.3, 1.6, 1.7], 
               [2, 2, 2, 2, 2]])
    
    

calc = Morisons(50, 5, 0, 0.1, 5)
z = calc.z_array(50, 0.1)
t = calc.time(10, 0.1)
Re = np.array((-10, 1e1, 1e3, 5e5, 1e10))
vel = calc.horizontal_velocity(1, t, z)
Re = np.array([[1, 2, 3], [4, 5, 6]])

