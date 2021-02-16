# Hydrodynamic Load Analysis using Morisons Equation
This Library allows users to model forces acting on a cylindrical bottom mounted offshore structure. It has the following utilities:
- Wave Field Meshing
- Wave Field Modelling
- Calculation of forces acting on the structure

The following conditions must be fulfilled to ensure accuracy of results:
- Unidirectional waves (Multi-Directional waves will be catered for in the upcoming patches)
- 2D wave field
- Deep Water Wave Conditions
- Small Body Assumption to ensure Morison's Equations hold true. 

Instructions for Use:
- Velocity, Pressure and Elevation functions within the Wave_Field Class use spatial grid, at a particular time value. The spatial grid needs to be created by using the meshgrid command by meshing the x and z arrays.  
- When using class Morisons, the calculations can be for multiple time steps, however, the value of x needs to be fixed. We assume that the diameter is small enough compared to the wavelength so much so that the force application remains the same across the surface. 

The second version of this library has following limitations:
- Doesn't allow broadcasting for the k, x, time and omega variables.
- Seed value for the kfromw function is set at 0. 
- The Wave_Field Class Functions, particularly, velocity and pressure, calculate these properties at a particular time. To obtain time marching solutions, loop needs to be used along with these functions. 
These limitations are imposed to keep the first release simple. A patch will be released in the next couple of months to fix this. 
 
