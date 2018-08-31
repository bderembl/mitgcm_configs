#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf

plt.ion()

flag_var = 0 #0: eta, 1: u

dir0 = '../run/mnc_test_0004/'
file1 = 'state.0000000000.t001.nc'

f1 = netcdf.netcdf_file(dir0 + file1,'r')

ti = f1.variables['T'][:].copy()
xc = f1.variables['X'][:].copy()
yc = f1.variables['Y'][:].copy()

si_t = len(ti)

nt1 = 0
nt2 = 70

if flag_var == 0:
  eta0 = f1.variables['Eta'][nt1,:,:].copy()
  eta1 = f1.variables['Eta'][nt2,:,:].copy()
elif flag_var == 1:
  u0 = f1.variables['U'][nt1,0,:,:-1].copy()
  u1 = f1.variables['V'][nt2,0,:,:-1].copy()


plt.figure()
for nt in range(0,si_t):
  if flag_var == 0:
    eta1 = f1.variables['Eta'][nt,:,:].copy()
    psi = eta1 - eta0
    vmax = 1e-3
  elif flag_var == 1:
    u1 = f1.variables['U'][nt,0,:,:-1].copy()
    psi = u1 - u0
    vmax = np.max(np.abs(psi))

  vmin = -vmax
  psi[0,0] = vmin
  psi[1,0] = vmax
  
  plt.contourf(xc,yc,psi,30,cmap=plt.cm.seismic,vmin=vmin,vmax=vmax)
  plt.draw()
  plt.colorbar()
  input("hit space");
  plt.clf()
