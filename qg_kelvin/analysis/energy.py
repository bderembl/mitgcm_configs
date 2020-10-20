#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit
import scipy.io.netcdf as netcdf
plt.ion()

flag_tile = 1
dir0 = '/run/media/bderembl/workd/MITgcm/myrun/test_kw_energetics/run04/'
#dir0 = '/home/bderembl/work/MITgcm/myrun/test_kw_energetics/run/'
dir1 = dir0 + 'mnc*/'
dir2 = dir0 + 'mnc_test_0001/'

file0 = 'grid.t*'
file1 = 'state.*'

if flag_tile == 0:
  file0 = 'grid.t001.nc'
  file1 = 'state.0000000000.t001.nc'
  

alphat = 2e-4
go = 9.81

# grid
if flag_tile == 1:
  f0 = mit.mnc_files(dir1 + file0)
else:
  f0 = netcdf.netcdf_file(dir2 + file0,'r')

RC    = f0.variables['RC'][:].copy()
DRC   = f0.variables['drC'][:].copy()
DRF   = f0.variables['drF'][:].copy()
RF    = f0.variables['RF'][:].copy()

XC    = f0.variables['XC'][:,:].copy()
YC    = f0.variables['YC'][:,:].copy()

si_y,si_x = XC.shape
si_z      = RC.size

dx = XC[1,1] - XC[0,0]
dy = YC[1,1] - YC[0,0]
dz = RC[1] - RC[0]

dv = np.abs(dx*dy*dz)

if flag_tile == 1:
  f1 = mit.mnc_files(dir1 + file1)
else:
  f1 = netcdf.netcdf_file(dir2 + file1,'r')

T = f1.variables['T'][:].copy()
si_t = len(T)

temp = f1.variables['Temp'][0,:si_z,:si_y,:si_x].copy()

RC3d = (0*temp + 1)*(RC.reshape(si_z,1,1))
# compute KE, PE
ener = np.zeros((si_t,5))

ener_tmp = np.zeros((si_t,si_z))
for nt in range (0,si_t):
  u    = f1.variables['U'][nt,:si_z,:si_y,:si_x].copy()
  v    = f1.variables['V'][nt,:si_z,:si_y,:si_x].copy()
#  w    = f1.variables['W'][nt,:si_z,:si_y,:si_x]
  temp = f1.variables['Temp'][nt,:si_z,:si_y,:si_x].copy()
  eta = f1.variables['Eta'][nt,:si_y,:si_x].copy()


  ener[nt,0] = 0.5*np.sum((u**2 + v**2 ))
  # without upper and lower boundaries
  ener[nt,1] = go*np.sum((1-alphat*(temp[1:-1,:,:]+273.1))*(RC[1:-1].reshape(si_z-2,1,1)))
  # with upper and lower boundaries
  #ener[nt,1] = go*np.sum((1-alphat*(temp[:,:,:]+273.1))*(RC[:].reshape(si_z,1,1)))
  #ener[nt,1] = go*np.sum((1-alphat*(temp[1:-1,:,:]+273.1))*(RC3d+eta)[1:-1,:,:])
  
#  ener_tmp[nt,:] = go*np.sum(np.sum((1-alphat*(temp+273.1))*(RC3d+eta),1),1)

#  ener[nt,1] = 0.5*go*np.sum((1-alphat*(temp+273.1))*(RF[:-1]**2-RF[1:]**2).reshape(si_z,1,1))



ener[:,0] = ener[:,0]*dv
ener[:,1] = ener[:,1]*dv
#ener[:,1] = ener[:,1]*dx*dy

enera = ener - ener[0,:]

plt.figure()
plt.plot(enera[:-1,0],'k',label='KEa')
plt.plot(enera[:-1,1],'r',label='PEa')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc=2)
plt.grid()
plt.xlabel('Time (day)')
plt.ylabel('Energy (m^5/s^2)')

#plt.savefig('energy_diag.pdf',bbox_inches='tight')
#np.savetxt('ener_mit.dat',enera)
