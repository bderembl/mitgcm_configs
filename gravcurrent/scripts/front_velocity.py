#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
from scipy.interpolate import interp1d

plt.ion()

dir0 = '../run/mnc_test_'
file1 = 'state.0000000000.t001.nc'

if len(sys.argv) > 1:
  dir0 = dir0 + str(format(sys.argv[1])).zfill(4) + '/'

# physical parameters
gg = 9.81       # m/s^2
sbeta = 7.4e-4  # psu^-1

f1 = netcdf.netcdf_file(dir0 + file1,'r')

tt = f1.variables['T'][:].copy()
xx = f1.variables['X'][:].copy()
zz = f1.variables['Z'][:].copy()

si_t = len(tt)
si_x = len(xx)
si_z = len(zz)

s = f1.variables['S'][:,:,:].copy().squeeze()
u = f1.variables['U'][:,:,:].copy().squeeze()

# parameters of the experiment
s1 = s[0,1,1]   # salty water
s0 = s[0,-1,-1] # fresh water

# binary image
s_bin = np.where(s<1e-3*s1,0.,1. )

# position of the front with the binary image
ipos = si_x - np.argmax(s_bin[:,-1,::-1],axis=1)
ipos = ipos[ipos<si_x]
imax = len(ipos)
xf = xx[ipos]

# parameters of the experiment
L0 = xx[ipos[0]]
H0 = -zz[-1] - zz[0]
gp = gg*sbeta*(s1-s0)
uref = np.sqrt(gp*H0)
tref = L0/uref

# compute velocity
dtu = np.diff(tt)
tu = 0.5*(tt[1:imax] + tt[:imax-1])
uf = np.diff(xf)/dtu[:imax-1]

uf2 = np.zeros(imax)
xf2 = np.zeros(imax)
nw = 10
for ii in range(0,imax):
  # maximum velocity near the front
  uf2[ii] = np.max(u[ii,-1,(ipos[ii]-nw):ipos[ii]])

  # find position of the front with an interpolation
  s_epsilon = np.linspace(1e-7,0,nw+1)
  fi = interp1d(s[ii,-1,(ipos[ii]-nw):ipos[ii]+1] + s_epsilon, xx[ipos[ii]-nw:ipos[ii]+1], kind='linear')
  xf2[ii] = fi(0.2*s1)

uf3 = np.diff(xf2)/dtu[:imax-1]

plt.figure()
#plt.loglog(tu,uf) 
plt.loglog(tt[1:imax]/tref,uf2[1:]/uref,'r',label='max velocity') 
plt.loglog(tu[1:imax]/tref,uf3[1:]/uref,'k',label='front velocity') 

plt.loglog(tu[2:]/tref,uf3[-1]/uref*(tu[2:]/tu[-1])**-(1/3),'k--',label='t^-1/3')
plt.xlabel('time')
plt.ylabel('velocity')
plt.legend()

