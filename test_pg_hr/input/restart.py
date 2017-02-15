#!/usr/bin/env python

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import pupynere as netcdf
import MITgcmutils as mit
from scipy import interpolate
from scipy.ndimage.filters import gaussian_filter

plt.ion()

binprec = '>f8'

flag_interp = 1


# dir0 = './'
# file1 = 'state.0000000000.t001.nc'

dir0 = './run/'
dir1 = dir0 + 'mnc_test*/'
file1 = 'state.*'

#f1 = netcdf.netcdf_file(dir0 + file1,'r')
f1 = mit.mnc_files(dir1 + file1)

T = f1.variables['T'][:]
nt = len(T)-1

uvel  = f1.variables['U'   ][nt,:,:,:]
vvel  = f1.variables['V'   ][nt,:,:,:]
theta = f1.variables['Temp'][nt,:,:,:]
eta   = f1.variables['Eta' ][nt,:,:]


si_z,si_y,si_x = theta.shape

uvel = uvel[:si_z,:si_y,:si_x]
vvel = vvel[:si_z,:si_y,:si_x]

if flag_interp:
  # old grid
  xx = np.linspace(0,1,si_x)
  yy = np.linspace(0,1,si_y)
  
  xog,yog = np.meshgrid(xx,yy)
  
  # new grid
  si_xn = 720
  si_yn = 720
  si_zn = si_z
  
  
  Lx = 3000.0e3
  Ly = 3000.0e3
  
  dx = Lx/si_xn;
  dy = Ly/si_yn;
  
  dx1 = dx*np.ones((si_xn))
  dy1 = dy*np.ones((si_yn))
  
  
  dx1.astype(binprec).tofile('dx.box')
  dy1.astype(binprec).tofile('dy.box')
  
  xn = np.linspace(0,1,si_xn)
  yn = np.linspace(0,1,si_yn)
  
  xng,yng = np.meshgrid(xn,yn)
  
  # fix theta on topo
  theta2 = 1.0*theta
  theta2[:,:,0]  = theta2[:,:,1]
  theta2[:,-1,:] = theta2[:,-2,:]
  
  
  uvel_n  = np.zeros((si_z,si_yn,si_xn))
  vvel_n  = np.zeros((si_z,si_yn,si_xn))
  theta_n = np.zeros((si_z,si_yn,si_xn))
  
  eta_n = np.zeros((si_yn,si_xn))
  
  for nz in range(0,si_z):
    fint = interpolate.interp2d(xx, yy,uvel[nz,:,:], kind='linear')
    uvel_n[nz,:,:] = fint(xn,yn)
  
    fint = interpolate.interp2d(xx, yy,vvel[nz,:,:], kind='linear')
    vvel_n[nz,:,:] = fint(xn,yn)
  
    fint = interpolate.interp2d(xx, yy,theta2[nz,:,:], kind='linear')
    theta_n[nz,:,:] = fint(xn,yn)
  
  fint = interpolate.interp2d(xx, yy,eta, kind='linear')
  eta_n = fint(xn,yn)
  
  # TODO: vertical interpolation
else:
  uvel_n  = uvel
  vvel_n  = vvel
  theta_n = theta
  eta_n   = theta
  

uvel_n.astype(binprec).tofile('uinit.box')
vvel_n.astype(binprec).tofile('vinit.box')
theta_n.astype(binprec).tofile('tinit.box')

eta_n.astype(binprec).tofile('einit.box')
