#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit
from scipy import interpolate
import sys

plt.ion()

# generate new grid in mygendata first!!

if len(sys.argv) > 1:
  iexp = int(sys.argv[1])
else:
  print("usage: python restart.py iexp")
  print("with iexp the number of the experiment")
  sys.exit(1)


binprec = '>f4'

dir0 = '../run/mnc_test_' + str(format(iexp)).zfill(4) + '/'
file1 = 'state.*'

#f1 = netcdf.netcdf_file(dir0 + file1,'r')
f1 = mit.mnc_files(dir0 + file1)

T = f1.variables['T'][:]
nt = len(T)-1

xp1 = f1.variables['Xp1'][:]
yp1 = f1.variables['Xp1'][:]
z = -f1.variables['Z'][:]
zl = -f1.variables['Zl'][:]

Lx = xp1[-1]
Ly = yp1[-1]
Lz = 2*z[-1]-zl[-1]

si_x = len(xp1) - 1
si_y = len(yp1) - 1
si_z = len(z)

uvel  = f1.variables['U'   ][nt,:,:,:]
vvel  = f1.variables['V'   ][nt,:,:,:]
theta = f1.variables['Temp'][nt,:,:,:]
eta   = f1.variables['Eta' ][nt,:,:]

# new grid
dxn = np.fromfile("dx.box",dtype=binprec, count=-1, sep='')
dyn = np.fromfile("dy.box",dtype=binprec, count=-1, sep='')
dzn = np.fromfile("dz.box",dtype=binprec, count=-1, sep='')
si_xn = len(dxn)
si_yn = len(dyn)
si_zn = len(dzn)

zc = np.cumsum(dzn)-0.5*dzn

zind1 = np.zeros(si_zn,dtype=int)
zind2 = np.zeros(si_zn,dtype=int)
wei1 = np.zeros(si_zn)
wei2 = np.zeros(si_zn)

for iz in range(0,si_zn):
  if zc[iz]<z[0]:
    zind1[iz] = 0
    zind2[iz] = 0
    wei1[iz] = 1
    wei2[iz] = 0
  elif zc[iz]>z[-1]:
    zind1[iz] = si_z-1
    zind2[iz] = si_z-1
    wei1[iz] = 0
    wei2[iz] = 1
  else:
    pos = np.where(z-zc[iz]<0,100000,z-zc[iz])
    zind2[iz] = np.argmin(np.abs(pos))
    zind1[iz] = zind2[iz] - 1
    wei1[iz] = -(z[zind2[iz]]-zc[iz])/(z[zind1[iz]]-z[zind2[iz]])
    wei2[iz] =  (z[zind1[iz]]-zc[iz])/(z[zind1[iz]]-z[zind2[iz]])

uvel = uvel[:si_z,:si_y,:si_x]
vvel = vvel[:si_z,:si_y,:si_x]

if si_xn != si_x:
  # old grid
  xx = np.linspace(0,1,si_x)
  yy = np.linspace(0,1,si_y)
  
  xog,yog = np.meshgrid(xx,yy)  
  
  #new grid
  xn = np.linspace(0,1,si_xn)
  yn = np.linspace(0,1,si_yn)
  
  xng,yng = np.meshgrid(xn,yn)
  
  # fix theta on walls
  theta2 = 1.0*theta
  theta2[:,:,0]  = theta2[:,:,1]
  theta2[:,-1,:] = theta2[:,-2,:]
  for nz in range(0,si_z):
    tz = np.mean(theta2[nz,:,:])
    theta2[nz,:,:] = np.where(theta2[nz,:,:] == 0,tz,theta2[nz,:,:])

  uvel_t  = np.zeros((si_z,si_yn,si_xn))
  vvel_t  = np.zeros((si_z,si_yn,si_xn))
  theta_t = np.zeros((si_z,si_yn,si_xn))
  
  eta_n = np.zeros((si_yn,si_xn))
  
  for nz in range(0,si_z):
    fint = interpolate.interp2d(xx, yy,uvel[nz,:,:], kind='linear')
    uvel_t[nz,:,:] = fint(xn,yn)
  
    fint = interpolate.interp2d(xx, yy,vvel[nz,:,:], kind='linear')
    vvel_t[nz,:,:] = fint(xn,yn)
  
    fint = interpolate.interp2d(xx, yy,theta2[nz,:,:], kind='linear')
    theta_t[nz,:,:] = fint(xn,yn)
  
  fint = interpolate.interp2d(xx, yy,eta, kind='linear')
  eta_n = fint(xn,yn)
  
  uvel_n  = np.zeros((si_zn,si_yn,si_xn))
  vvel_n  = np.zeros((si_zn,si_yn,si_xn))
  theta_n = np.zeros((si_zn,si_yn,si_xn))

  for nz in range(0,si_zn):
    uvel_n[nz,:,:] = wei1[nz]*uvel_t[zind1[nz],:,:] + wei2[nz]*uvel_t[zind2[nz],:,:]
    vvel_n[nz,:,:] = wei1[nz]*vvel_t[zind1[nz],:,:] + wei2[nz]*vvel_t[zind2[nz],:,:]
    theta_n[nz,:,:] = wei1[nz]*theta_t[zind1[nz],:,:] + wei2[nz]*theta_t[zind2[nz],:,:]

else:
  uvel_n  = uvel
  vvel_n  = vvel
  theta_n = theta
  eta_n   = eta
  

uvel_n.astype(binprec).tofile('uinit.box')
vvel_n.astype(binprec).tofile('vinit.box')
theta_n.astype(binprec).tofile('tinit.box')
eta_n.astype(binprec).tofile('einit.box')
