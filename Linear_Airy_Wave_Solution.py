C

"""
Created on Wed Feb  10 22:22:10 2021

@author: goharshoukat
"""
import scipy.optimize
import numpy as np
import math
import sys
import matplotlib.pyplot as plt

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
        return (scipy.optimize.fsolve(solve_for_k, 0))

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


#Chile Class of the Set_Wave_Field Parent. 
#Pressure, Velocity fields will be calculated and plotted here
#Linear Airy Solution

class Wave_Field(Set_Wave_Field):
    
    def __init__(self, depth, time_period, tide_velocity, amplitude, \
                 x: np.ndarray, delta_z):
        super().__init__(depth, time_period, tide_velocity, amplitude)
        self.x = x
        self.z = np.arange(-depth, 0 + delta_z, delta_z)
        self.X, self.Z = np.meshgrid(self.x, self.z)

    def meshgrid(self):
        return (self.X, self.Z)
    
    #wave surface eleveation function. 
    #two dimensional time varying array
    #row = z variation
    #column = time variation
    #x position is a float value, not an array
    #code can deal with arrays of K and omega, but works with only one array 
    #either x or time. user will have to define which one to use per run. 
    
    #dimensions of k if array is passed = m x 1
    #dimension of x if array is passed = 1 x m
    #dimensions of omega if array is passed = m x 1
    #dimension of t if array is passed = 1 x m
    
    def elevation(self, x, time):
        if type(x) == np.ndarray and type(time) == np.ndarray:
            raise TypeError ('x and time can not both be arrays')
        else:
            return   (self.amplitude * np.cos(self.kfromw()\
                                                * x - self.omega*time))    
                
    def velocity_potential_spatial_variable(self, time: float):
        #for evolution with time, run a for loop with different values of Time
        #The argument need to be fixed in the next patch. Function needs to 
        #be able to either use omega from self or as an argument passed directly
        #on to the function
        k = self.kfromw()
        return(self.amplitude * self.gravity / self.omega * \
            np.cosh(k * (self.Z + self.depth))/np.cosh(k * self.depth) * \
                np.sin(k * self.X -self.omega * time))
        
    def horizontal_velocity_spatial_variable(self, time: float):
        k = self.kfromw()
        return (self.amplitude * self.omega * \
                (np.cosh(k * (self.Z + self.depth)))/np.sinh(k * self.depth) \
                    *np.cos(k*self.X - self.omega * time))
    
    def vertical_velocity_spatial_variable(self, time: float):
        k = self.kfromw()        
        return (self.amplitude * self.omega * \
                (np.sinh(k * (self.Z + self.depth)))/np.sinh(k * self.depth) \
                    *np.sin(k*self.X - self.omega * time))
        
    def linearPressure_Dynamic(self, time: float, rho):
        
        k = self.kfromw()
        return (self.amplitude * rho * self.gravity * \
                (np.cosh(k * (self.Z + self.depth)))/np.sinh(k * self.depth) \
                    *np.cos(k*self.X - self.omega * time))
    
    #we pass both the horizontal and vertical velocity. the function plots the 
    #correct angle and mangnitude of the quivers. 
    def plot(self, pressure, hor_velocity, vert_velocity, time):
        plt.figure()
        plt.contourf(self.X, self.Z, pressure)
        plt.colorbar()
        plt.quiver(self.X, self.Z, \
                   self.horizontal_velocity_spatial_variable(time), \
                       self.vertical_velocity_spatial_variable(time))
        plt.plot(self.x, self.elevation(self.x, time))
        plt.xlabel('x [m]')
        plt.ylabel('z [m]')
        plt.title('T = ' + str(time) + 's')
        plt.show()