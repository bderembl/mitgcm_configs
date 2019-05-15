#!/usr/bin/env python

import glob,os,re
import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
from scipy import interpolate
import MITgcmutils as mit

plt.ion()

binprec = '>f4'

flag_interp = 1

#dir0 = '../run/mnc_test_0002/'
dir0 = '../run/run0003/'

file1 = 'state.*'

file_u = 'U*'
file_v = 'V*'
file_t = 'T*'
file_s = 'S*'
file_e = 'Eta*'

flag_nc = 0;

if flag_nc:
  allfiles = sorted(glob.glob(dir0 + file1));
  si_t  = len(allfiles);
  
  f1 = netcdf.netcdf_file(allfiles[-1],'r')
  
  T = f1.variables['T'][:].copy()
  nt = len(T)-1
  
  uvel  = f1.variables['U'   ][nt,:,:,:].copy()
  vvel  = f1.variables['V'   ][nt,:,:,:].copy()
  theta = f1.variables['Temp'][nt,:,:,:].copy()
  salt  = f1.variables['S'   ][nt,:,:,:].copy()
  eta   = f1.variables['Eta' ][nt,:,:].copy()
  
  xx   = f1.variables['X'][:].copy()
  yy   = f1.variables['Y'][:].copy()

  si_z,si_y,si_x = theta.shape
  
  uvel = uvel[:si_z,:si_y,:si_x]
  vvel = vvel[:si_z,:si_y,:si_x]

  Lx = xx[-1] + xx[0]
  Ly = yy[-1] + yy[0]

else:
  iters = mit.mds.scanforfiles(dir0 + file_u)
  si_t = len(iters)
  
  # grid files
  xx = mit.rdmds(dir0 + 'XC').squeeze()
  yy = mit.rdmds(dir0 + 'YC').squeeze()
  
  
  nt = -1
  
  uvel  = mit.rdmds(dir0 + file_u, itrs = [iters[nt]]).squeeze()
  vvel  = mit.rdmds(dir0 + file_v, itrs = [iters[nt]]).squeeze()
  theta = mit.rdmds(dir0 + file_t, itrs = [iters[nt]]).squeeze()
  salt  = mit.rdmds(dir0 + file_s, itrs = [iters[nt]]).squeeze()
  eta   = mit.rdmds(dir0 + file_e, itrs = [iters[nt]]).squeeze()

  si_z,si_y,si_x = theta.shape

  Lx = xx[0,-1] + xx[0,0]
  Ly = yy[-1,0] + yy[0,0]


if flag_interp:
  # old grid
  xx = np.linspace(0,1,si_x)
  yy = np.linspace(0,1,si_y)
  
  xog,yog = np.meshgrid(xx,yy)
  
  # new grid
  si_xn = 800
  si_yn = 800
  si_zn = si_z
    
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

  salt2 = 1.0*salt
  salt2[:,:,0]  = salt2[:,:,1]
  salt2[:,-1,:] = salt2[:,-2,:]
  
  
  uvel_n  = np.zeros((si_z,si_yn,si_xn))
  vvel_n  = np.zeros((si_z,si_yn,si_xn))
  theta_n = np.zeros((si_z,si_yn,si_xn))
  salt_n  = np.zeros((si_z,si_yn,si_xn))
  
  eta_n = np.zeros((si_yn,si_xn))
  
  for nz in range(0,si_z):
    fint = interpolate.interp2d(xx, yy,uvel[nz,:,:], kind='linear')
    uvel_n[nz,:,:] = fint(xn,yn)
  
    fint = interpolate.interp2d(xx, yy,vvel[nz,:,:], kind='linear')
    vvel_n[nz,:,:] = fint(xn,yn)
  
    fint = interpolate.interp2d(xx, yy,theta2[nz,:,:], kind='linear')
    theta_n[nz,:,:] = fint(xn,yn)
  
    fint = interpolate.interp2d(xx, yy,salt2[nz,:,:], kind='linear')
    salt_n[nz,:,:] = fint(xn,yn)
  
  fint = interpolate.interp2d(xx, yy,eta, kind='linear')
  eta_n = fint(xn,yn)
  
  # TODO: vertical interpolation
else:
  uvel_n  = uvel
  vvel_n  = vvel
  theta_n = theta
  salt_n  = salt
  eta_n   = eta
  

uvel_n.astype(binprec).tofile('uinit.box')
vvel_n.astype(binprec).tofile('vinit.box')
theta_n.astype(binprec).tofile('tinit.box')
salt_n.astype(binprec).tofile('sinit.box')

eta_n.astype(binprec).tofile('einit.box')
