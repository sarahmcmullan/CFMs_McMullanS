################################# Catastrophic Fragmentation Models Input ##################################
#
# Created by Sarah McMullan, s.mcmullan16@imperial.ac.uk
#
# Last updated: 26/07/2017
#
############################################################################################################
#
### Description:
#
# Input file for calculating three different continuous fragmentation models: 
#		- Pancake Model, Chyba et al., 1993, (PM)
#		- Debris Cloud Model, Hills and Goda, 1993, (DC) 
#		- Chain Reaction Model, Avramenko et al., 2014, (CR)
#
# There are various options that can be implemented in the models:
#		- Atmospheric Density Model:
#			- 'Isothermal'
#			- 'CurveFit': Curve-fit equation for the 1976 U.S. Standard Atmosphere,from Wheeler et al., 2017 
#			- 'Interpolation': Interpolation from the 1976 U.S. Standard Atmosphere, requires data every 10 m
#		- Increasing Strenth with decreasing fragment size (only applies to Pancake and Debris Cloud):
#			- 'Yes': active (no longer require f_p)
#			- 'No': inactive, strength remains constant
#		- Mass loss due to ablation affecting radius (only applies to Pancake and Debris Cloud):
# 			- 'Yes': active
#			- 'No': inactive
# 		- Energy deposited:
#			- 'Released': total energy released
#			- 'Deposited': total energy released minus the energy radiated
#		- Coefficients used for each model:
# 			- 'Combined': coefficients defined used for all three models are the same
#			- 'Individual': coefficients are indiviually defined for all three models      
# If a coefficent is not in use for a specific run, either comment it out or set it to 'NA'
#
# ??????Assumptions made:
# 	- Earth can be assumed to be flat, valid if theta isn't really shallow
#	- Gravity and lift are negligable, valid for all medium to large meteoroids 
#
# Code is vaild when:
#   - theta >~ 15 degrees
#	- radius of meteoroid >~ 5 m 
#	- max altitude <= 86km
#???????
# Symbols used:
# 	v = meteorite speed (m s^-1)
#	z = altitude (m)
#	rho_m = density of meteorite (kg m^-3)
# 	theta = atmospheric entry angle (degrees)
# 	rad = radius of meteoroid (m^2)
# 	rho_a_model = select which atmospheris density model to implement
# 	t = time (secs)
# 	mstrength = tensile strength of meteoroid (Pa), ~1e6-1e7 Pa for stony meteoroids 
# 	C_D = drag coefficient
#	C_H = ablation parameter (kg J^-1)
#	Q = heat of ablation
#	C_L = leift coefficent
# 	C_F = constant, related to fragmentation
# 	C_R = proportion of emergy released from ablation
#
######################################################################################################

print 'Continuous Fragmentation Models Start'

########################################## Import libraries ##########################################

import numpy as n
import matplotlib.pyplot as plt
import sys
import os
from CFModels_Earth import *

######################################################################################################

######################################## Select Model to Use  ########################################

m = CFModels( )

######################################################################################################

######################################### Editable parameters ########################################

### Meteorite initial conditions
m.v0 = 19040. 		# initial meteorite speed immediately before entering atmosphere (m s^-1)
m.z0 = 86000.		# starting altitude (m), must be <= 86 km
m.rho_m = 3300.		# density of meteoroid (kg m^-3)
theta0_deg = 17.		# inital angle of trajectory (degrees)
m.theta0 = theta0_deg * (n.pi/180.) # convert trajectory to radians for calculations
m.rad0 = 9.95		# initial radius (m) 
m.rho_a_model = 'Interpolate'	# select which atmospheric density model to use: 'Isothermal', 'StandardEquation', 'Interpolate'

# Time parameters of model
m.t0 = 0.			# start time (secs)
m.tend = 20.		# end time (secs)
m.dt = 2.e-4		# timestep size

### Coefficients
m.combined_individual = 'Individual'  # Select whether to use 'Combined' or 'Individual' coefficients

### Combined Coefficients
m.C_D = 1.5
m.abp = 1.2e-8
m.C_L = 1.e-3
m.mstrength0 = 5.e5	# tensile strength of meteoroid (Pa)
m.energy_calc = 'Deposited'

### Chain Reaction Coefficients
m.CRC_D = 1.5	# drag coeff.
m.CRabp = 1.2e-8
m.CRC_L = 1.e-3
m.alpha = 0.177
m.CRmstrength0 = 4.e5
m.C_F = 1.3		# fragmentation coeff.
m.C_R = 0.37 #0.093
m.CRenergy_calc = 'Deposited'     # Selet energy calculation to use: 'Simple', 'Released', 'Deposited'

### Data file names
m.datadirectory = "Data_CR"		# location of where to sabe output data
m.dedzdata = "dEdZ"	    # name of output data file for dE/dz
m.runname = "Test"			# name of specific run
		

### U.S. Standard Atmosphere Data Import 
US_Stand_At_Loc = "1976AtmosphericDensityPlusH.csv"
data = n.loadtxt(US_Stand_At_Loc,delimiter=" ",skiprows=1)
m.zdata = data [:,0]
m.rho_a_data = data [:,1]
m.H_interpolator = data [:,2]

######################################################################################################

############################################# Run Models ##############################################

print 'Running Chain Reaction Model...'
m.RunCR( m.t0, m.tend)

######################################################################################################

############################################ Output data #############################################

print 'Saving Data'
m.SaveDEDz_Single('CR')

######################################################################################################

print 'Continuous Fragmentation Models End'

######################################################################################################


 
