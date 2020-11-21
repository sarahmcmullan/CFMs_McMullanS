import matplotlib.pyplot as plt
import numpy as np
import CFM_plot as p
import os


model_name = ['PM','DC','CR']
model_colour = ['#9900cc','#0099cc','#ff8000']
for i in range(len(model_name)):
  [z,E,t,v,m,r,theta] = np.loadtxt('Data/dEdZ{}v19040_r9_rho_m3300_theta17_sigma500000_abp1.200000e-08_CD1.5_CS0.3Interpolate.csv'.format(model_name[i]),unpack=True,skiprows = 1)#, usecols= (0,1,2,3,4,5))
  plt.figure(1)
  dEdZ = p.dEdZ_vs_z(E,z,'Yes',model_colour[i])
  plt.figure(2)
  spreading_ratio = r/r[0]
  rvz = p.x_vs_z(spreading_ratio,z,model_colour[i], 'Spreading Ratio', xmax = 12) #### Change to spreading ratio
  plt.figure(3)
  rvz = p.t_vs_x(t,spreading_ratio,model_colour[i], 'Spreading Ratio', ymax = 12) #### Change to spreading ratio


plt.tight_layout()
if not os.path.exists('Figures'):
  os.makedirs('Figures')
plt.figure(1)
plt.savefig('Figures/dEdz_v_z.pdf')
plt.figure(2)
plt.savefig('Figures/r_v_z.pdf')
plt.figure(3)
plt.savefig('Figures/t_v_r.pdf')

plt.show()

