import matplotlib.pyplot as plt
import numpy as np


def dEdZ_vs_z(dEdZ_data, altitude_data, plot_Chelyabinsk_Brown, plot_colour, xmax='not defined', ymax = 70, xmin = 0,ymin = 0):
  if plot_Chelyabinsk_Brown == 'Yes':
    data = np.loadtxt("/home/sm1216/asgard_bkp/sm1216/pancake_model/ElliotsCode/ElliotsCode/Output/Brown.dat",delimiter=",",skiprows=1)
    zbr = np.array(list(reversed(data[:,0])))		# create array of altutude from Brown et al.
    Ebr = np.array(list(reversed(data[:,1]))) 		# create array of energy deposition from Brown et al.
    z = np.linspace(0.,50.,401)
    Eb = np.interp(z,zbr,Ebr,left=0.,right=0.)		# interpolate Brown energy deposition to calulated 
    plt.plot(Eb,z,'k-',alpha=0.8,label="Brown et al. (2013)",lw=4)		# plot Brown et al. data
  
  plt.plot(dEdZ_data, altitude_data, color = plot_colour,lw=2)
  plt.grid(lw=0.5)
  plt.ylabel('Altitude (km)',fontsize=12)
  plt.xlabel('Energy Deposition (kT/km)',fontsize=12)
  if xmax == 'not defined':
    xmax = (np.max(dEdZ_data))*1.1
  plt.gca().set_xlim(left=xmin, right=xmax)
  plt.gca().set_ylim(bottom=ymin, top=ymax)

def x_vs_z(x_data, altitude_data, plot_colour, xlabel='...', xmax='not defined', ymax = 70, xmin = 0,ymin = 0):
  plt.plot(x_data, altitude_data, color = plot_colour,lw=2)
  plt.grid(lw=0.5)
  plt.ylabel('Altitude (km)',fontsize=12)
  plt.xlabel('{}'.format(xlabel),fontsize=12)
  if xmax == 'not defined':
    xmax = (np.max(x_data))*1.1
  plt.gca().set_xlim(left=xmin, right=xmax)
  plt.gca().set_ylim(bottom=ymin, top=ymax)

def t_vs_x(time_data, y_data,plot_colour, ylabel='...', xmax='not defined', ymax = 'not defined', xmin = 0,ymin = 0):
  plt.plot(time_data, y_data, color = plot_colour,lw=2)
  plt.grid(lw=0.5)
  plt.xlabel('Time (secs)',fontsize=12)
  plt.ylabel('{}'.format(ylabel),fontsize=12)
  if ymax == 'not defined':
    ymax = (np.max(y_data))*1.1
  if xmax == 'not defined':
    xmax = (np.max(time_data))
  plt.gca().set_xlim(left=xmin, right=xmax)
  plt.gca().set_ylim(bottom=ymin, top=ymax)


