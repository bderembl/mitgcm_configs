#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
import f90nml

plt.ion()

dir0 = '../run/'
file1 = 'state.0000000000.t001.nc'

if len(sys.argv) > 1:
  dir1 = dir0 + 'mnc_test_' + str(format(sys.argv[1])).zfill(4) + '/'

# physical parameters
gg = 9.81       # m/s^2
sbeta = 7.4e-4  # psu^-1

nml = f90nml.read(dir0+'data')
nmldiag = f90nml.read(dir0+'data.diagnostics')

dt = nml['parm03']['dumpfreq']


f1 = netcdf.netcdf_file(dir1 + file1,'r')

tt = f1.variables['T'][:].copy()
xx = f1.variables['X'][:].copy()
zz = f1.variables['Z'][:].copy()

si_t = len(tt)
si_x = len(xx)
si_z = len(zz)

s = f1.variables['S'][:,:,:,:].copy().squeeze()
u = f1.variables['U'][:,:,:,:-1].copy().squeeze()
w = f1.variables['W'][:,:,:,:].copy().squeeze()

# parameters of the experiment
s1 = s[0,1,1]   # salty water
s0 = s[0,-1,-1] # fresh water
ipos = si_x - np.argmax(s[0,-1,::-1])
L0 = xx[ipos]
H0 = -zz[-1] - zz[0]
Lt = xx[-1] - xx[0]
dx = xx[1] - xx[0] # assume uniform grid
dz = zz[0] - zz[1] # assume uniform grid
gp = gg*sbeta*(s1-s0)
uref = np.sqrt(gp*H0)
tref = L0/uref

# change z = 0 at the bottom
zz = H0 + zz


b = -gg*sbeta*s

ke = 0.5*(u**2 + w**2)
pe = -b*zz[np.newaxis,:,np.newaxis]

# vorticity at vertices
dudz = -np.diff(u,axis=1)[:,:,1:]/dz
dwdx = np.diff(w,axis=2)[:,1:,:]/dx

omega = dwdx - dudz

ke_s = ke.sum((1,2))*dx*dz
pe_s = pe.sum((1,2))*dx*dz
omega_s = omega.sum((1,2))*dx*dz

# plt.figure()
# plt.plot(tt/tref,ke_s/(uref**2*H0*L0),'k',label='ke')
# plt.plot(tt/tref,pe_s/(uref**2*H0*L0),'r',label='pe')
# plt.plot(tt/tref,(ke_s+pe_s)/(uref**2*H0*L0),'g--',label='tot')
# plt.xlabel('time')
# plt.ylabel('E')
# plt.legend()

plt.figure()
plt.loglog(tt/tref,ke_s/(uref**2*H0*L0),'k',label='ke')
plt.loglog(tt/tref,pe_s/(uref**2*H0*L0),'r',label='pe')
plt.loglog(tt/tref,(ke_s+pe_s)/(uref**2*H0*L0),'g--',label='tot')
plt.xlabel('time')
plt.ylabel('E')
plt.legend()
