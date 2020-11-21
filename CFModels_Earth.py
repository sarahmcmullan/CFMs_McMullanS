### Continuous Fragmentation Models
# Sarah McMullan 
#
# Description: Functions to run models described in:
#     - Chyba, 1993
#     - Hills & Goda, 1993
#     - Avramenko, 2014
# with edited equations to make models run as similar as possible apart from the change in radius equation
# - Multiple options for Atmospheric Model to use to run model with
#
################################# Continuous Fragmentation Models ####################################
#
# Created by Sarah McMullan, November 2016
# Last Edited 28/11/2018
#
######################################################################################################
#
import numpy as n
from scipy.integrate import ode 
import os

class CFModels( object ):
  # Defining model paramenters to be input
  v0 = None
  z0 = None
  rho_m = None
  theta0 = None
  mstrength0 = None
  f_p = None
  rad0 = None
  rho_a_model = None
  t0 = None
  tend = None
  dt = None
  combined_individual = None
  C_D = None
  abp = None
  C_L = None
  energy_calc = None
  PMC_D = None		
  PMC_H = None
  PMC_L = None
  PMmstrength0 = None
  PMenergy_calc = None
  PMincr_strength = None
  PMablation_rad = None
  DCC_D = None			
  DCC_S = None
  DCC_L = None
  DCmstrength0 = None
  DCenergy_calc = None
  DCincr_strength = None
  DCablation_rad = None
  CRC_D = None		
  CRC_H = None		
  C_F = None			
  C_R = None
  alpha = None
  CRC_L = None
  CRmstrength0 = None
  CRenergy_calc = None
  fragmentation_begun = 0

  # Earth Parameters
  H = 8000.
  g = 9.81
  rho_g = 1.386
  RE = 6378000.

  # Unit conversions
  j2kt = 2.39e-13
  deg2rad = n.pi/180.
  m2km = 1e-3

  # Create list in to which to store the raw data
  raw_results_PM = []  
  raw_results_DC = []  
  raw_results_CR = []  

  # Import data for interpolating Standard Atmospheric Density
  data = None 
  zdata = None 
  rho_a_data = None 
  H_interpolator = None 

  # Creates Class
  def __init__( self ): 
    ''''''

  # function to calculate atmospheric density with altitude
  def AtmosphericDensity( self, ht):
    if self.rho_a_model == 'Isothermal': 
      return self.rho_g * n.exp(-(ht/self.H))
    elif self.rho_a_model == 'CurveFit':
      return -140.2*n.exp(-0.000187*ht) + 141.4*n.exp(-0.000186*ht)
    elif self.rho_a_model == 'Interpolate':
      if ht % 10 != 0:
        zref = int(ht/10) # calculate reference number of point closest but smaller than z
        rho_a = (self.rho_a_data[zref]*n.exp(-(ht-self.zdata[zref])/self.H_interpolator[zref]))
      else:
        rho_a = (self.rho_a_data[int(ht/10)])
      return rho_a     

  # function to calculate the cross-sectional area assuiming spherical
  def AreaRadius( self, rad0):
    A0 = n.pi * rad0**2
    return A0

  # function to calculate the mass assuiming spherical
  def MassCalc( self, rad0):
    m0 = self.rho_m * (4./3.) * n.pi * (self.rad0)**3
    return m0

  # function of all Pancake Model ODEs
  def PM_ODEs( self, t, W):
    rad = W[0]; m = W[1]; v = W[2]; z = W[3]; y = W[4]; theta = W[5]; E = W[6];

    if self.combined_individual == 'Combined':
      C_D = self.C_D
      abp = self.abp
      C_L = self.C_L
      mstrength0 = self.mstrength0
      energy_calc = self.energy_calc
    else:
      C_D = self.PMC_D
      abp = self.PMabp
      C_L = self.PMC_L
      mstrength0 = self.PMmstrength0
      energy_calc = self.PMenergy_calc

    A = n.pi * rad**2

    fdrag = -(C_D * self.AtmosphericDensity(z) * A * v**2)/2.

    dv = fdrag/m + self.g*n.sin(theta)

    dm = -(abp * self.AtmosphericDensity(z) * A * v**3)/2

    dtheta = ((self.g * n.cos(theta))/v) - ((self.C_L * self.AtmosphericDensity(z) * A * v)/(2 * m)) - ((v * n.cos(theta))/(self.RE + z))

    dz = -v * n.sin(theta)

    Pram = (self.AtmosphericDensity(z)* v**2)

    if self.PMincr_strength == 'Yes':
      mfr = (9. * n.pi * m**3.)/(16. * A**3. * self.rho_m**2.)
      m0 = self.MassCalc(self.rad0)
      mstrength = mstrength0 * (mfr/m0)**-self.alpha
    else:
      mstrength = mstrength0 

    if self.PMablation_rad == 'Yes': 			# does ablation effect radius
      if Pram < mstrength:
        dy = 0.
        drad = ((1./3.) * (rad/m) * dm)
      else:
        dy = ((C_D*self.AtmosphericDensity(z)*(v**2))/(2*rad*self.rho_m))
        drad = y + ((1./3.) * (rad/m) * dm)		

    else:
      if Pram < mstrength:
        dy = 0.
      else:
        dy = ((C_D*self.AtmosphericDensity(z)*(v**2))/(2*rad*self.rho_m))	

      if rad > (self.f_p * self.rad0):
        drad = 0.
      else:
        drad = y

    if self.PMenergy_calc == 'Released':
      dE = (C_D * A *((self.AtmosphericDensity(z) * n.abs(v)**3)/2.)) + (abp * A *((self.AtmosphericDensity(z) * n.abs(v)**5)/(4.)))
    elif self.PMenergy_calc == 'Deposited':
      dE = (C_D * A *((self.AtmosphericDensity(z) * n.abs(v)**3.)/2.)) + ((1. - self.C_R)*(abp * A *((self.AtmosphericDensity(z) * n.abs(v)**5.)/(4.))))
    else:
      dE = 0.
    return [ drad, dm, dv, dz, dy, dtheta, dE]



  # function of all Debri Cloud Model ODEs
  def DC_ODEs( self, t, W):
    rad = W[0]; m = W[1]; v = W[2]; z = W[3]; theta = W[4]; E = W[5]; 

    global fragmentation_begun

    if self.combined_individual == 'Combined':
      C_D = self.C_D
      abp = self.abp
      C_L = self.C_L
      mstrength0 = self.mstrength0
    else:
      C_D = self.DCC_D
      abp = self.DCabp
      C_L = self.DCC_L
      mstrength0 = self.DCmstrength0
  
    A = n.pi * rad**2

    fdrag = -(C_D * self.AtmosphericDensity(z) * A * v**2)/2.

    dv = fdrag/m + self.g*n.sin(theta)

    dm = -(abp * self.AtmosphericDensity(z) * A * v**3)/2

    dtheta = ((self.g * n.cos(theta))/v) - ((self.C_L * self.AtmosphericDensity(z) * A * v)/(2 * m)) - ((v * n.cos(theta))/(self.RE + z))

    dz = -v * n.sin(theta)

    Pram = (self.AtmosphericDensity(z)* v**2)
    if t == 0.:
       fragmentation_begun = 0

    # Increasing Strength
    if self.DCincr_strength == 'Yes':
      mfr = (9. * n.pi * m**3.)/(16. * A**3. * self.rho_m**2.)
      m0 = self.MassCalc(self.rad0)
      mstrength = mstrength0 * (mfr/m0)**-self.alpha
    else:
      mstrength = mstrength0

    Vmax = (mstrength/1.293)**(1./2.)

    if self.DCablation_rad == 'Yes':				# does ablation affect radius
      if Pram < mstrength:
        drad = 0. + ((1./3.) * (rad/m) * dm)
      else:
        drad = (((7./2.)*self.C_S*(self.AtmosphericDensity(z)/self.rho_m))**(1./2.) * v) + ((1./3.) * (rad/m) * dm)
    else:
      if Pram < mstrength and fragmentation_begun == 0:
        drad = 0.
      elif Pram >= mstrength:
        drad = ((7./2.)*self.DCC_S*(self.AtmosphericDensity(z)/self.rho_m))**(1./2.) * v
        fragmentation_begun = 1
      elif Pram < mstrength and self.fragmentation_begun == 1:
        drad = ((7./2.)*self.DCC_S*(self.AtmosphericDensity(z)/self.rho_m))**(1./2.) * v
        if v > Vmax:
          fragmentation_begun = 1
        else:
          fragmentation_begun = 2
      else:
        drad = 0.

    if self.DCenergy_calc == 'Released':
      dE = (C_D * A *((self.AtmosphericDensity(z) * n.abs(v)**3)/2.)) + (abp * A *((self.AtmosphericDensity(z) * n.abs(v)**5)/(4.)))
    elif self.DCenergy_calc == 'Deposited':
      dE = (C_D * A *((self.AtmosphericDensity(z) * n.abs(v)**3.)/2.)) + ((1. - self.C_R)*(abp * A *((self.AtmosphericDensity(z) * n.abs(v)**5.)/(4.))))
    else:
      dE = 0.
    
    return [ drad, dm, dv, dz, dtheta, dE]


  # function of all Chain Reaction Model ODEs
  def CR_ODEs( self, t, W):
    m = W[0]; v = W[1]; z = W[2]; A = W[3]; E = W[4]; theta = W[5];
	
    if self.combined_individual == 'Combined':
      C_D = self.C_D
      abp = self.abp
      C_L = self.C_L
      mstrength0 = self.mstrength0
    else:
      C_D = self.CRC_D
      abp = self.CRabp
      C_L = self.CRC_L
      mstrength0 = self.CRmstrength0

    fdrag = -(C_D * self.AtmosphericDensity(z) * A * v**2)/2.

    dv = fdrag/m + self.g*n.sin(theta)

    dm = -(abp * self.AtmosphericDensity(z) * A * v**3)/(2)

    dtheta = ((self.g * n.cos(theta))/v) - ((self.C_L * self.AtmosphericDensity(z) * A * v)/(2 * m)) - ((v * n.cos(theta))/(self.RE + z))

    dz = -v * n.sin(theta)

    Pram = (self.AtmosphericDensity(z)* v**2)

    if self.CRenergy_calc == 'Released':
      dE = (C_D * A *((self.AtmosphericDensity(z) * n.abs(v)**3)/2.)) + (abp * A *((self.AtmosphericDensity(z) * n.abs(v)**5)/(4.)))
    elif self.CRenergy_calc == 'Deposited':
      dE = (C_D * A *((self.AtmosphericDensity(z) * n.abs(v)**3.)/2.)) + ((1. - self.C_R)*(abp * A *((self.AtmosphericDensity(z) * n.abs(v)**5.)/(4.))))
    else:
      dE = 0.

    mfr = ((9. * n.pi * m**3.)/(16. * A**3. * self.rho_m**2.))
    m0 = self.MassCalc(self.rad0)
    mstrength = mstrength0 * (mfr/m0)**-self.alpha

    if Pram < mstrength:
      dA = (2./3.) * (A/m) * dm

    else:
      Afr = self.C_F * (((Pram - mstrength)**(1./2.))/((m**(1./3.))*(self.rho_m**(1./6.)))) * A    
      dA = ((2./3.) * (A/m) * dm) + Afr

    return [ dm, dv, dz, dA, dE, dtheta]

  # run model solving Pancake Model ODEs for dt
  def RunPM( self, t0, tend):
    # Calculate Mass
    m0 = self.MassCalc(self.rad0)
    # set up differentials
    W0 = [ self.rad0, m0, self.v0, self.z0, 0., self.theta0, 0.]
    r = ode(self.PM_ODEs).set_integrator('dopri5')#,etol=1e-6)
    r.set_initial_value( W0, t0).set_f_params()
    # calculate total energy of system
    etot = 0.5 * (m0 * (self.v0**2))
    # set up lists to store parameter values			
    vl = [W0[2]]
    zl = [W0[3]]
    ml = [W0[1]]
    thetal = [W0[5]]
    Al = [n.pi*(W0[0]**2)]
    radl = [W0[0]]
    el = [W0[6]]
    dedz = []
    time = []
    # do ODE calulations
    while r.successful() and r.t < tend and r.y[3] >= 0:
      time += [r.t]
      #print 'r time = ', r.t
      r.integrate(r.t+self.dt)
      vl += [r.y[2]]
      zl += [r.y[3]]
      ml += [r.y[1]]
      thetal += [r.y[5]]
      if self.PMenergy_calc == 'Simple':
        el += [0.5 * ((m0 * self.v0**2) - (r.y[1] * r.y[2]**2))]
      else:        
        el += [r.y[6]]
      dedz += [(el[-1] - el[-2]) / (zl[-1] - zl[-2])]
      radl += [r.y[0]]
      Al += [n.pi * (r.y[0]**2)]
    # save all lists to single list of arrays
    self.raw_results_PM = [n.array(time), n.array(vl), n.array(zl), n.array(ml), n.array(el), n.array(dedz), n.array(radl), n.array(Al), n.array(etot), n.array(thetal)]


  # run model solving Debris Cloud Model ODEs for dt
  def RunDC( self, t0, tend):
    # Calculate Mass
    m0 = self.MassCalc(self.rad0)
    # set up differentials
    W0 = [ self.rad0, m0, self.v0, self.z0, self.theta0, 0.]
    r = ode(self.DC_ODEs).set_integrator('dopri5')
    r.set_initial_value( W0, t0).set_f_params()
    # calculate total energy of system
    etot = 0.5 * (m0 * (self.v0**2))
    # set up lists to store parameter values			
    vl = [W0[2]]
    zl = [W0[3]]
    ml = [W0[1]]
    thetal = [W0[4]]
    Al = [n.pi*(W0[0])**2.]
    radl = [W0[0]]
    el = [W0[5]]
    dedz = []
    time = []
    # Say that fragmentation hasn't begun at the beginning of run
    fragmentation_begun = 0
    # do ODE calulations
    while r.successful() and r.t < tend and r.y[3] >= 0:
      time += [r.t]
      r.integrate(r.t+self.dt)
      vl += [r.y[2]]
      zl += [r.y[3]]
      ml += [r.y[1]]
      thetal += [r.y[4]]
      if self.DCenergy_calc == 'Simple':
        el += [0.5 * ((m0 * self.v0**2) - (r.y[1] * r.y[2]**2))]
      else:        
        el += [r.y[5]]
      dedz += [(el[-1] - el[-2]) / (zl[-1] - zl[-2])]
      radl += [r.y[0]]
      Al += [n.pi*(r.y[0])**2.]
    # save all lists to single list of arrays
    self.raw_results_DC = [n.array(time), n.array(vl), n.array(zl), n.array(ml), n.array(el), n.array(dedz), n.array(radl), n.array(Al), n.array(etot), n.array(thetal)]#, n.array(Praml), n.array(mstrengthl)]

  # run model solving Chain Reaction Model ODEs for dt
  def RunCR( self, t0, tend):
    # Calculate Mass
    m0 = self.MassCalc(self.rad0)
    A0 = self.AreaRadius( self.rad0)

    # set up differentials
    W0 = [ m0, self.v0, self.z0, A0, 0., self.theta0]
    r = ode(self.CR_ODEs).set_integrator('dopri5')
    r.set_initial_value( W0, t0).set_f_params()
    # calculate total energy of system
    etot = 0.5 * (m0 * (self.v0**2))
    # set up lists to store parameter values			
    vl = [W0[1]]
    zl = [W0[2]]
    ml = [W0[0]]
    Al = [W0[3]]
    thetal = [W0[5]]
    radl = [self.rad0]
    el = [W0[4]]
    dedz = []
    time = []
    # do ODE calulations
    while r.successful() and r.t < tend and r.y[3] >= 0:
      time += [r.t]
      r.integrate(r.t+self.dt)
      vl += [r.y[1]]
      zl += [r.y[2]]
      ml += [r.y[0]]
      if self.CRenergy_calc == 'Simple':
        el += [0.5 * ((m0 * self.v0**2) - (r.y[0] * r.y[1]**2))]
      else:        
        el += [r.y[4]]
      thetal += [r.y[5]]
      dedz += [(el[-1] - el[-2]) / (zl[-1] - zl[-2])]
      Al += [r.y[3]]
      radl += [(r.y[3]/n.pi)**(1./2.)]
    # save all lists to single list of arrays
    self.raw_results_CR = [n.array(time), n.array(vl), n.array(zl), n.array(ml), n.array(el), n.array(dedz), n.array(radl), n.array(Al), n.array(etot), n.array(thetal)]


  # Output time array
  def Time_All(self):
    return [self.raw_results_PM[0],self.raw_results_DC[0],self.raw_results_CR[0]]
	
  # Output velocity array
  def Velocity_All(self):
    return [self.raw_results_PM[1],self.raw_results_DC[1],self.raw_results_CR[1]]

  # Output altitude array
  def Altitude_All(self):
    return [self.raw_results_PM[2],self.raw_results_DC[2],self.raw_results_CR[2]]

  # Output mass array
  def Mass_All(self):
    return [self.raw_results_PM[3],self.raw_results_DC[3],self.raw_results_CR[3]]
	
  # Output energy release array	
  def EnergyRelease_All(self):
    return [self.raw_results_PM[4],self.raw_results_DC[4],self.raw_results_CR[4]]

  # Output de/dz array
  def DeDz_All(self):
    return [self.raw_results_PM[5],self.raw_results_DC[5],self.raw_results_CR[5]]

  # Output radius array
  def MeteoroidRadius_All(self):
    return [self.raw_results_PM[6],self.raw_results_DC[6],self.raw_results_CR[6]]

  # Output area array
  def MeteoroidArea_All(self):
    return [self.raw_results_PM[7],self.raw_results_DC[7],self.raw_results_CR[7]]

  # Output total energy released value
  def TotalEnergy_All(self):
    return [self.raw_results_PM[8],self.raw_results_DC[8],self.raw_results_CR[8]]

  # Output trajectory array
  def Trajectory_All(self):
    return [self.raw_results_PM[9],self.raw_results_DC[9],self.raw_results_CR[9]]

  # Output time array for an individual model
  def Time_Single(self,model_name):
    if model_name == 'PM':
      return self.raw_results_PM[0]
    elif model_name == 'DC':
      return self.raw_results_DC[0]
    elif model_name == 'CR':
      return self.raw_results_CR[0]
	
  # Output velocity array for an individual model
  def Velocity_Single(self,model_name):
    if model_name == 'PM':
      return self.raw_results_PM[1]
    elif model_name == 'DC':
      return self.raw_results_DC[1]
    elif model_name == 'CR':
      return self.raw_results_CR[1]

  # Output altitude array for an individual model
  def Altitude_Single(self,model_name):
    if model_name == 'PM':
      return self.raw_results_PM[2]
    elif model_name == 'DC':
      return self.raw_results_DC[2]
    elif model_name == 'CR':
      return self.raw_results_CR[2]

  # Output mass array for an individual model
  def Mass_Single(self,model_name):
    if model_name == 'PM':
      return self.raw_results_PM[3]
    elif model_name == 'DC':
      return self.raw_results_DC[3]
    elif model_name == 'CR':
      return self.raw_results_CR[3]
	
  # Output energy release array for an individual model	
  def EnergyRelease_Single(self,model_name):
    if model_name == 'PM':
      return self.raw_results_PM[4]
    elif model_name == 'DC':
      return self.raw_results_DC[4]
    elif model_name == 'CR':
      return self.raw_results_CR[4]

  # Output de/dz array for an individual model
  def DeDz_Single(self,model_name):
    if model_name == 'PM':
      return self.raw_results_PM[5]
    elif model_name == 'DC':
      return self.raw_results_DC[5]
    elif model_name == 'CR':
      return self.raw_results_CR[5]

  # Output radius array for an individual model
  def MeteoroidRadius_Single(self,model_name):
    if model_name == 'PM':
      return self.raw_results_PM[6]
    elif model_name == 'DC':
      return self.raw_results_DC[6]
    elif model_name == 'CR':
      return self.raw_results_CR[6]

  # Output area array for an individual model
  def MeteoroidArea_Single(self,model_name):
    if model_name == 'PM':
      return self.raw_results_PM[7]
    elif model_name == 'DC':
      return self.raw_results_DC[7]
    elif model_name == 'CR':
      return self.raw_results_CR[7]

  # Output total energy released value for an individual model
  def TotalEnergy_Single(self,model_name):
    if model_name == 'PM':
      return self.raw_results_PM[8]
    elif model_name == 'DC':
      return self.raw_results_DC[8]
    elif model_name == 'CR':
      return self.raw_results_CR[8]

  # Output trajectory array for an individual model
  def Trajectory_Single(self,model_name):
    if model_name == 'PM':
      return self.raw_results_PM[9]
    elif model_name == 'DC':
      return self.raw_results_DC[9]
    elif model_name == 'CR':
      return self.raw_results_CR[9]


  # Save the output data arrays for all three models
  def SaveDEDz_All(self):
    if not os.path.exists(self.datadirectory):
      os.makedirs(self.datadirectory)

    modelname = ['PM','DC','CR']
    for i in range(len(modelname)):
      runname = "{}v{}_r{}_rho_m{}_theta{}_sigma{}_abp{:.6e}_CD{}_CS{}".format(modelname[i],int(self.v0),int(self.rad0),int(self.rho_m),int(self.theta0/self.deg2rad),int(self.mstrength0),self.abp,round(self.C_D,6),round(self.DCC_S,6))
      n.savetxt('{}/{}{}{}.csv'.format(self.datadirectory,self.dedzdata,runname,self.rho_a_model), n.c_[self.Altitude_All()[i][:-1:10]*self.m2km, -self.DeDz_All()[i][::10]*(self.j2kt/self.m2km), self.Time_All()[i][::10], self.Velocity_All()[i][:-1:10], self.Mass_All()[i][:-1:10], self.MeteoroidRadius_All()[i][:-1:10], self.Trajectory_All()[i][:-1:10]/self.deg2rad], fmt='%.18e', delimiter = " ", newline='\n', header='Altitude    Energy Deposition per unit height   Time    Velocity    Mass    Radius    Trajectory', footer='', comments='# ')


  # Save the output data arrays for individual models
  def SaveDEDz_Single(self, model_name):
    if not os.path.exists(self.datadirectory):
      os.makedirs(self.datadirectory)

    if model_name == 'PM':
      runname = "PMv{}_r{}_rho_m{}_theta{}_sigma{}_abp{:.6e}_CD{}_CS{}".format(int(self.v0),int(self.rad0),int(self.rho_m),int(self.theta0/self.deg2rad),int(self.mstrength0),self.abp,round(self.C_D,6),round(self.DCC_S,6))
      n.savetxt('{}/{}{}{}.csv'.format(self.datadirectory,self.dedzdata,runname,self.rho_a_model), n.c_[self.Altitude_Single('PM')[:-1:10]*self.m2km, -self.DeDz_Single('PM')[::10]*(self.j2kt/self.m2km), self.Time_Single('PM')[::10], self.Velocity_Single('PM')[:-1:10], self.Mass_Single('PM')[:-1:10], self.MeteoroidRadius_Single('PM')[:-1:10], self.Trajectory_Single('PM')[:-1:10]/self.deg2rad], fmt='%.18e', delimiter = " ", newline='\n', header='Altitude    Energy Deposition per unit height   Time    Velocity    Mass    Radius    Trajectory', footer='', comments='# ')
    elif model_name == 'DC':
      runname = "DCv{}_r{}_rho_m{}_theta{}_sigma{}_abp{:.6e}_CD{}_CS{}".format(int(self.v0),int(self.rad0),int(self.rho_m),int(self.theta0/self.deg2rad),int(self.mstrength0),self.abp,round(self.C_D,6),round(self.DCC_S,6))
      n.savetxt('{}/{}{}{}.csv'.format(self.datadirectory,self.dedzdata,runname,self.rho_a_model), n.c_[self.Altitude_Single('DC')[:-1:10]*self.m2km, -self.DeDz_Single('DC')[::10]*(self.j2kt/self.m2km), self.Time_Single('DC')[::10], self.Velocity_Single('DC')[:-1:10], self.Mass_Single('DC')[:-1:10], self.MeteoroidRadius_Single('DC')[:-1:10], self.Trajectory_Single('DC')[:-1:10]/self.deg2rad], fmt='%.18e', delimiter = " ", newline='\n', header='Altitude    Energy Deposition per unit height   Time    Velocity    Mass    Radius    Trajectory', footer='', comments='# ')
    elif model_name == 'CR':
      runname = "CRv{}_r{}_rho_m{}_theta{}_sigma{}_abp{:.6e}_CD{}_CS{}".format(int(self.v0),int(self.rad0),int(self.rho_m),int(self.theta0/self.deg2rad),int(self.mstrength0),self.abp,round(self.C_D,6),round(self.DCC_S,6))
      n.savetxt('{}/{}{}{}.csv'.format(self.datadirectory,self.dedzdata,runname,self.rho_a_model), n.c_[self.Altitude_Single('CR')[:-1:10]*self.m2km, -self.DeDz_Single('CR')[::10]*(self.j2kt/self.m2km), self.Time_Single('CR')[::10], self.Velocity_Single('CR')[:-1:10], self.Mass_Single('CR')[:-1:10], self.MeteoroidRadius_Single('CR')[:-1:10], self.Trajectory_Single('CR')[:-1:10]/self.deg2rad], fmt='%.18e', delimiter = " ", newline='\n', header='Altitude    Energy Deposition per unit height   Time    Velocity    Mass    Radius    Trajectory', footer='', comments='# ')








































