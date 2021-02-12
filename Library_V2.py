#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  10 22:22:10 2021

@author: goharshoukat
"""
import scipy.optimize
import numpy as np
import math

class Set_Wave_Field:
    #All parameters passed should be in SI Units
    gravity = 9.81
    
    #takes the basic parameters of wave field
    def __init__(self,depth, time_period, tide_velocity, amplitude):
        self.depth = depth
        self.tide_velocity = tide_velocity
        self.omega = 2*math.pi/time_period
        self.amplitude = amplitude
        
    #determines the wave number iteratively by solving the dispersion relation
    def kfromw(self):
        def solve_for_k(k):
            return (self.tide_velocity * k + np.sqrt(self.gravity * k * \
                                        np.tanh(k * self.depth)) - self.omega)     
        return (scipy.optimize.fsolve(solve_for_k,0))

    #calculates the wave length
    def wave_length(self):        
        return (2 * math.pi/self.kfromw())
    
    #determines water depth regime. Equations change for deep, shallow and 
    #intermediate
    def depth_regime(self):
        ratio = self.depth / self.wave_length()
        if ratio <= 1/20:
            print ("Shallow Water")
     
        elif ratio >= 1/2:
            print ("Deep Water")

        elif ratio > 1/20 and ratio < 1/2:
            print ("Intermediate Water")
            
    #redundancy check. depth_regime and mu determine the water depth regime
    #mu > O(1) for linear Airy solution to work
    def mu(self):
        ratio = self.kfromw()*self.depth
        if ratio > 1:
            print("Deep Water")
        else:
            print('Water not Deep Enough')
        return None
        
    #Determines if Morison's Equation is applicable or not. 
    #Makes use of the displacement parameter for judgement. 
    #KC can also be used. KC = 2pi A/D
    def diffraction_effects(self, diameter):
        ratio = diameter/self.wave_length()
        if ratio < 0.6:
            print('Small Body Assumption')
        else:
            print('Large Body Assumption')
        return None
    
    #Calculation of wave steepness parameter to limit wave field
    #within linear airy wave theory
    def epsilon(self):
        return (self.kfromw() * self.amplitude)
        
    