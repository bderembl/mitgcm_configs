#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit
import scipy.io.netcdf as netcdf
import sys

plt.ion()

iexp = ['6']

dir0 = '../run10_p'
dir1 = dir0 + iexp[0] + '/mnc_test_*/'
#file_s = 'state*.nc'
file_m = 'mean*.nc'
file_g = 'grid*.nc'

#f1 = netcdf.netcdf_file(dir1 + file_m,'r')

fg = mit.mnc_files(dir1 + file_g)
hfacc = fg.variables['HFacC'][:,:,:]
hfacw = fg.variables['HFacW'][:,:,:]
hfacs = fg.variables['HFacS'][:,:,:]

xx = fg.variables['X'][:]
yy = fg.variables['Y'][:]
xp1 = fg.variables['Xp1'][:]
yp1 = fg.variables['Yp1'][:]
zz = -fg.variables['Z'][:]
zl = -fg.variables['Zl'][:]

fg.close()

xg,yg = np.meshgrid(xx,yy) 

Lx = xp1[-1]
Ly = yp1[-1]
Lz = 2*zz[-1]-zl[-1]

si_x = len(xx)
si_y = len(yy)
si_z = len(zz)

dx = xx[1] - xx[0]
dy = yy[1] - yy[0]

uvel_me  = np.zeros((si_z,si_y,si_x))
vvel_me  = np.zeros((si_z,si_y,si_x))
theta_me = np.zeros((si_z,si_y,si_x))

nme = 0

for  ie in range(0,len(iexp)):

  dir1 = dir0 + iexp[ie] + '/mnc_test_*/'
  f1 = mit.mnc_files(dir1 + file_m)
  
  it = f1.variables['iter'][:]
  si_t = len(it)
  
  
  it0 = 0
  it1 = si_t-1
  for it in range(it0,it1):
  
    uvel  = f1.variables['UVEL' ][it,:,:,:]
    vvel  = f1.variables['VVEL' ][it,:,:,:]
    theta = f1.variables['THETA'][it,:,:,:]
  
    uvel_me  += uvel[:si_z,:si_y,:si_x]
    vvel_me  += vvel[:si_z,:si_y,:si_x]
    theta_me += theta[:si_z,:si_y,:si_x]
    nme += 1
  
  f1.close()
  
uvel_me  /= nme
vvel_me  /= nme
theta_me /= nme
  
# create netcdf file
fileout = dir0 + iexp[0] + '/average.nc'
f = netcdf.netcdf_file(fileout,'w')

f.createDimension('iter',1)
f.createDimension('Z',si_z)
f.createDimension('Y',si_y)
f.createDimension('X',si_x)

ito = f.createVariable('iter', 'f', ('iter',))
zpo = f.createVariable('Z', 'f', ('Z',))
ypo = f.createVariable('Y', 'f', ('Y',))
xpo = f.createVariable('X', 'f', ('X',))

uo  = f.createVariable('U' , 'f', ('iter','Z','Y','X',))
vo  = f.createVariable('V' , 'f', ('iter','Z','Y','X',))
to  = f.createVariable('Temp' , 'f', ('iter','Z','Y','X',))

ito[0] = 1
zpo[:] = zz
ypo[:] = yy
xpo[:] = xx


uo[:,:,:] = uvel_me[:,:,:]
vo[:,:,:] = vvel_me[:,:,:]
to[:,:,:] = theta_me[:,:,:]


f.close()
print ("nb points in file: {0}".format(nme))
