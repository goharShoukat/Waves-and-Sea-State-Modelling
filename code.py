# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 12:13:37 2020

@author: Gohar Shoukat
\
    
"""
import math
import numpy as np
import API_RP_2A_WSD_2000_Ed21_S3_2007 as lib
import matplotlib.pyplot as plt
#import pandas as pd

#Data inpub
#All data is in SI Units
#Design significant wave height 
Hs = 10
#H=2 times the amplitude
H= 1.8*Hs
#density of water
rho = 1024
#Design Peak=absolute period
Tp=14
#cero-rossing period
Tz=11
#pile diameter
dia=6
#water depth
h = 20
#total tower length
tower_length=90
#gravity 
g=9.81
#design current speed - colineal current
u_c = 1
kinematic_visc = 1.787e-6
#phi - angle between the waves and current
phi = 0
#tide
tide = 4

w_a = lib.absolute_frequency(Tp)
#solve relative frequency and max horizontal velocity due to wave only
#marine growth 0-46m as per API standard
growth = 3.81e-2
dia_adjusted = dia + 2*growth
#k is the diffraction parameter,
#L is the wave length

k, L = lib.diffraction_parameter(g, h, tide, u_c, w_a)
#w_r=math.sqrt(g*k*math.tanh(k*h))
#T_r=2*math.pi/w_r


#Call Diffraction Function to see if Structure is large or small
#and if diffraction plays an important role or not
#False for small structure
boolean = lib.diffraction_theory(dia, L)

#time for which simulation needs to be run in sec
start_time = 0
total_time = 14
#delta_time, difference between time steps
delta_time = 0.1
#define time matrix for solving equations for a set duration of the time
#the second argument has delta_time added to account for the function 
#limitation. np.arange terminates before the stop value. 
time = np.arange(start_time,total_time+delta_time,delta_time)
#size of time function for further operations is needed. 
length_time = len(time)

#determine the depth regime
ratio, x = lib.depth(L, h, None)

#stepping value for the height
delta_z = 0.1



u0, z, s, sigma, accel0 = \
lib.horizontal_velocity_acceleration_amplitude(g, H, Tp, k, phi, tide, x, L, h, 
                                               delta_z, H/2,w_a, time)
                                             
                                                                                      
#Final velocity at each step after addition with current
#size of array for further initializations
'''
plt.plot(u0, z)
# naming the x axis 
plt.xlabel('Amplitude of Velocity') 
# naming the y axis 
plt.ylabel('z') 
'''
#velocity to be calculated at intervals of delta
#size function used to determine the size of array
#z = 0 @ surface and z = max at seabed
#np.where(z>s,z,0)


#create a two dimensional array for depth and time. depth on y axis and 
#time on x axis. nested loop for solving it. 
#use the length functions to initialize the array
#each column in the u_horizontal_velocity,u_hor corresponds with 
#a specif depth. likewise, each row is a time step. 
#u_hor is the horizontal velocity only
#u is the velocity after the impact of current is added

acceleration = np.zeros((len(z),length_time))
u = np.zeros((len(z),length_time))
for j in range(length_time):
    for i in range(len(z)):
        #compare values of U against the displacement equation and remove 
                                             #the values where there isnt any 
        if z[i]<=s[j]:
            u[i,j] = u_c + u0[i,j] * math.cos(-w_a*time[j])
            acceleration[i,j] = accel0[i,j] * math.sin(-w_a*time[j])                                             #water
        else:

            u[i,j] = 0
            acceleration[i,j] = 0

#use the above caclculations for determining the Constants
#Re changes throughout the depth. We will update it at each iteration. 
#Re = np.zeros(length_u, length_time)

Re = (u_c+u0)*dia_adjusted/kinematic_visc
C_d = C_d2 = np.zeros((len(z),length_time))
CM = CM2 = np.zeros((len(z),length_time))
#Relative Roughness e
e = growth / dia_adjusted
#############
#############
#############
'''
values for the sea surface
'''
#############
#############
######check############
#change this to Tp
KC_surface = u0[z==0,:]*Tp/dia_adjusted
r_surface = u_c/(u0[z==0,:]+u_c)
theta_star_surface = math.pi/2-np.arctan(-1*r_surface/np.sqrt(1-r_surface))
correction_factor_surface = (1+r_surface)*2*theta_star_surface/math.pi
KC_corrected_surface = KC_surface*correction_factor_surface

#read from figure using e from page 64
CDS = 1.05
#KCr is corrected KC as per page 65
KCr_CDS = KC_corrected_surface/CDS
#read Cd/CDS from the graph
CD_CDS_surface = 1.5
#drag coefficient
C_d[z==0] = CD_CDS_surface*CDS


#for CM, do the following
#from page 66, read the value of CM against KCr/CDS
CM[z==0] = 1.4
#############
#############
#############


#############
#############
#############
'''
values at the bottom
'''
KC_seabed = u0[z==(-h-tide),:]*Tp/dia_adjusted
r_seabed = u_c/(u0[z==(-h-tide),:]+u_c)
theta_star_seabed = math.pi/2-np.arctan(-r_seabed/np.sqrt(1-r_seabed**2))
correction_factor_seabed = (1+r_seabed)*2*theta_star_seabed/math.pi
KC_corrected_seabed = KC_seabed*correction_factor_seabed

  
#KCr is corrected KC as per page 65
KCr_CDS_seabed = KC_corrected_seabed/CDS
#read Cd/CDS from the graph
CD_CDS_seabed = 1.6 
#drag coefficient
C_d[z==(-h-tide)] = CD_CDS_seabed*CDS

#for CM, do the following
#from page 66, read the value of CM against KCr/CDS
CM[z==(-h-tide)] = 1.5 

#replace all C_d values above z=0 with c_d @ z=0

#############
#############
#############
#Now, interpolate the Cd and CM over the depth of the sea

#interpolate between top and bottom using depth. 
#find gradient
m_Cd = (C_d[z==(-h-tide)]-C_d[z==0])/(-h-tide)
m_CM = (CM[z==(-h-tide)] - CM[z==0])/(-h-tide)
#finding KC at each value using equation of line


#interpolating C_d and CM between surface and seabed only. 
#following identifies the index for seabed and surface and runs only \
#between these two points. 

indx_surface = np.where(z==0)[0].tolist()
indx_seabed = np.where(z==(-h-tide))[0].tolist()
for i in range(indx_seabed[0],indx_surface[0]):
    C_d2[i] = m_Cd*z[i] + C_d[z==0]
    CM2[i] = m_CM*z[i] + CM[z==0]
    
for i in range(indx_surface[0],len(z)):
    C_d2[i] = C_d2[z==0]
    CM[i] = CM[z==0]

#C_d2 = np.where(z<0,C_d,C_d[z==0,:])
#CM2 = np.where(z<0,CM,CM[z==0])


#force calculation at each time step at each depth
    #shear force
    #recheck the magnitude of u in the equation, something is wrong here 
    #possibly

force = np.zeros((len(z), length_time))
force_drag = np.zeros((len(z), length_time))
force_M = np.zeros((len(z), length_time))

for j in range(length_time):
    for i in range(len(z)):
        if z[i]<=s[j]:
   
            force_drag[i,j] =  0.5 * rho * dia_adjusted\
            *abs(u[i,j])*(u[i,j]) * C_d2[i,j] 
        
            force_M [i,j] = rho * CM2[i,j] * math.pi\
            * (dia_adjusted**2) / 4 * acceleration[i,j]
                
# 0.5 * rho * dia_adjusted *abs(u[i,j])*(u[i,j]) * C_d2[i]             
#total force.
force = force_drag + force_M
from scipy import integrate
#
#Integration along the depth at each time step
#integrate each column. 
I = np.zeros((length_time))
I_force_Cd = np.zeros((length_time))
I_force_M = np.zeros((length_time))
for i in range(length_time):
    I[i] = np.trapz(force[:,i],z)
    I_force_Cd[i] = integrate.simps(force_drag[:,i],z)
    I_force_M[i] = integrate.simps(force_M[:,i],z)

#part c of qs 1
max_drag = max(I_force_Cd)
min_drag = min(I_force_Cd)
max_inertia = max(I_force_M)
min_inertia = min(I_force_M)
max_total_force = max(I)
min_total_force = min(I)
plt.plot(time, I)
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(time,I, label = "Total Force")
ax.plot(time,I_force_Cd, label = 'Drag Force')
ax.plot(time,I_force_M, label = 'Inertial Force')
ax.legend()
plt.xlabel('Time [s]')
plt.ylabel('Force [N/m]')
plt.title('Total Force')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.savefig('G:\REM\Semester 2\Environmental conditions for marine\
            renewable\Total Force Simps.png')
plt.close()
#Moment Diagrams
moment = np.zeros((len(z),length_time))
for j in range(length_time):
    for i in range(len(z)):
        moment [i,j] = force[i,j] * (z[i]+h+tide) 

I_moment = np.zeros(length_time)
for i in range(length_time):
    I_moment[i] = integrate.simps(moment[:,i],z)

max(I_moment)
min(I_moment)
plt.plot(time, I_moment)
plt.xlabel('Time [s]')
plt.ylabel('Moment [N/m . m]')
plt.title('Moment')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.savefig('G:\REM\Semester 2\Environmental conditions for\
            marine renewable\Moment.png')
plt.close()


##########################################################################
##########################################################################
#Using image processing instead of linear interpolation for graphical data
##########################################################################
##########################################################################


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#                                Part B
#######################################################################
#######################################################################
#######################################################################
########################################################################
#wave celerity
c = 1.092 * g * Tp / ( 2 * math.pi)
#c=14.5297
#(breaking waves)
beta = 0.48 
#beta = 0.70 (broken waves)
#u_impact is the impact velocity
u_impact = beta * c

#impact pressure
#breaking waves 
k1 = 5.98
p = k1 * u_impact**2 * rho

#inline force/unit length
fmax = 2 * math.pi * rho * dia_adjusted/2 * c**2

#impact duration
Ti = 0.406 * dia_adjusted /(2 * c)

#curling factor
lambda_ = 0.5
nu = 0.95 * H

time2 = np.arange(0,1.5*Ti,0.0001)
tprime = time2 - dia_adjusted/(c*64)

#function values
f = np.zeros(len(time2))

for i in range(len(time2)):
    if time2[i] < (dia_adjusted/2/c/8) and time2[i]>=0:
       f[i]= rho*dia_adjusted/2*(c**2)*\
       ((2*math.pi-(2*np.sqrt(c/dia_adjusted/2*time2[i])*\
                    np.arctanh(np.sqrt(1-(c/dia_adjusted/2*time2[i]*0.25))))))
       
    elif tprime[i] >= 0.046875*dia_adjusted/c and tprime[i] <= 0.1875*dia_adjusted/c:
         f[i]=rho*dia_adjusted/2*c*c*\
         (math.pi*np.sqrt(dia_adjusted/2/6/c/tprime[i])-\
          (8*c*tprime[i]/3/dia_adjusted/2)**0.25*\
          np.arctanh(np.sqrt(1-c*tprime[i]/dia_adjusted/2*\
                             np.sqrt(6*c*tprime[i]/dia_adjusted/2))))

    else:
        f[i] = 0
        
f[0] = fmax
F = f * lambda_ * nu

Moment_breaking_waves = F * (h+tide+H/2)
plt.plot(time2,F)
plt.plot(time2,Moment_breaking_waves)

#Grid Convergence Analysis. 
#Export data of dz = 0.01 first
'''
df = pd.DataFrame({'time':time,'I':I})
filepath = 'G:\REM\Semester 2\Environmental conditions fo\
r marine renewable\dz'+ str(delta_z) + '.xlsx' 
df.to_excel(filepath, index = False)'''
#define time matrix for solving equations for a set duration of the time
#the second argument has delta_time added to account for the function 
#limitation. np.arange terminates before the stop value. 
#import data for dz=0.1

import pandas as pd
df = pd.read_excel(r'G:\REM\Semester 2\Environmental conditions\
                   for marine renewable\dz0.1.xlsx')
plt.plot(df['time'],df['I'])

plt.plot(time,I)
max(I) - max(df['I'])

fig_gridIndepence = plt.figure()
ax = plt.subplot(111)
ax.plot(time,I, label = "dz = dt = 1")
ax.plot(df['time'],df['I'], label = 'dz = dt = 0.1')
ax.legend()
plt.xlabel('Time [s]')
plt.ylabel('Force [N/m]')
plt.title('Total Force against time with different grid sizes')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.savefig('G:\REM\Semester 2\Environmental conditions\
            for marine renewable\Grid Independence.png')
plt.close()

'''time3 = np.arange(start_time,total_time+delta_time,delta_time)
force = data[:,1]
plt.plot(time3,force)
plt.plot(time, I)'''