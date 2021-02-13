# Hydrodynamic Load Analysis using Morisons Equation
This Library allows users to model forces acting on a cylindrical bottom mounted offshore structure. It has the following utilities:
- Wave Field Meshing
- Wave Field Modelling
- Calculation of forces acting on the structure

The following conditions must be fulfilled to ensure accuracy of results:
- Unidirectional waves (Multi-Directional waves will be catered for in the upcoming patches)
- Deep Water Wave Conditions
- Small Body Assumption to ensure Morison's Equations hold true. 

The second version of this library has following limitations:
- Doesn't allow broadcasting for the k, x, time and omega variables.
- Seed value for the kfromw function is set at 0. 

These limitations are imposed to keep the first release simple. A patch will be released in the next couple of months to fix this. 
 
