"""
A freely-propagating, premixed flat flame with multicomponent
transport properties.
"""

###############################################################################
#import 
###############################################################################
import cantera as ct
import numpy as np
import csv 
###############################################################################
#programme
###############################################################################
gas = ct.Solution('LiDryer.cti')

# Set initial states
p = 101325  # pressure [Pa]
Tin = 300.0  # temperature [K]
reactants = 'H2:2, O2:1, N2:3.76'  
gas.TPX = Tin, p, reactants

# Simulation parameters
width = 0.015  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

# Set up flame object
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=2, slope=0.006, curve=0.006, prune=0.00001)
f.set_grid_min(1e-8)
f.set_max_grid_points(1,10000)

#############################################
# Solve with mixture-averaged transport model
f.transport_model = 'Mix'
f.soret_enabled = True
f.solve(loglevel=loglevel, auto=True)
print("Mix: ", "p =",p, "T =", Tin,  "sL =", f.u[0], "sb =",f.u[1403])

# write the velocity, temperature, density, and mole fractions to a CSV file
###############################################################################

###############
# plot results
###############
#Plot the velocity, temperature, density
z = f.flame.grid
T = f.T
u = f.u

f.write_csv('phi=1.0.csv',species='Y',quiet=False)