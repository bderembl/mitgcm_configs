#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf

plt.ion()


dir0 = '../run/mnc_test_0019/'
file1 = 'state.0000000000.t001.nc'
file2 = 'grid.t001.nc'
file3 = 'phiHydLow.0000000000.t001.nc'

f1 = netcdf.netcdf_file(dir0 + file1,'r')
f2 = netcdf.netcdf_file(dir0 + file2,'r')
f3 = netcdf.netcdf_file(dir0 + file3,'r')

topo = f2.variables['Depth'][:].copy().squeeze()

zl = f1.variables['Zl'][:].copy().squeeze()
xx = f1.variables['X'][:].copy().squeeze()
t1 = f1.variables['T'][:].copy().squeeze()
t3 = f3.variables['T'][:].copy().squeeze()

si_x = len(xx)
si_z = len(zl)
si_t1 = len(t1)
si_t3 = len(t3)

L0 = 1000
u0z = 1e-3     # vertical shear (du/dz) s-1
U0 = u0z*L0
rho0 = 1e3
gg = 9.81      # gravity
alphaK = 2.0e-4
#H0 = 99.5 # m not exactly 100m to avoid discretization issue 

# mnc004
H0 = 49.5 # m not exactly 100m to avoid discretization issue 
rho0z = 9.2e-4   # stratification (d rho /dz)
nu = 1.0  # viscosity

# # mnc010
# H0 = 99.5 # m not exactly 100m to avoid discretization issue 
# rho0z = 1.1e-5   # stratification (d rho /dz)
# nu = 0.1 # viscosity

# mnc018
H0 = 99.5 # m not exactly 100m to avoid discretization issue 
rho0z = 1.1e-4   # stratification (d rho /dz)
nu = 1.0 # viscosity


J = gg*rho0z/(rho0*u0z**2)
S = H0/L0
delta = (nu/u0z/L0**2)**(1./3)
adimfac = delta*np.sqrt(J)*S**2/2

zl = (zl - zl[-1])/L0
xx = (xx - (xx[0] + xx[-1])/2)/L0
topo = (topo[0] - topo)/L0

xg, zg = np.meshgrid(xx,zl)

# adim time
t3 = t3*u0z

nt = -1

u = f1.variables['U'][nt,:,:].copy().squeeze()
v = f1.variables['V'][nt,:,:].copy().squeeze()
w = f1.variables['W'][nt,:,:].copy().squeeze()
temp = f1.variables['Temp'][nt,:,:].copy().squeeze()

u = u[:,:-1] # remove last point

#background profile
u0 = u[:,0]
ua = 0.*u 
for ix in range(0,si_x):
  ua[:,ix] = u[:,ix] - u0

# surface pressure
phiHydLow = f3.variables['phiHydLow'][:,:].copy().squeeze()
Pb = phiHydLow/U0**2

nx0 = 50
nx1 = 250


# no division by dx because we integrate just after
dhdx = np.diff(topo)
dPbdx = np.diff(Pb)

surf_pres_drag = np.sum((Pb[:,nx0+1:nx1]+Pb[:,nx0:nx1-1])/2*dhdx[nx0:nx1-1],1)
#surf_pres_drag = np.sum((Pb[:,1:]+Pb[:,:-1])/2*dhdx,1)
#surf_pres_drag = np.sum((Pb[:,1:])*dhdx*np.diff(xx),1)
#surf_pres_drag = np.sum(topo[1:]*dPbdx*np.diff(xx),1)

wavedrag = np.zeros((si_t1,si_z))

for nt in range(0,si_t1):
  wavedrag[nt,:] =  sum(ua[:,nx0+1:nx1].T*w[:,nx0+1:nx1].T*np.diff(xx[nx0:nx1]).reshape(nx1-nx0-1,1),0)



plt.figure()
plt.plot(t3,surf_pres_drag/adimfac)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('time')
plt.ylabel('Surface pressure drag')


plt.figure()
plt.plot(wavedrag[-1,:]/adimfac,zl,)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.xlim([-0.35,0.0])
plt.xlabel('wave drag')
plt.ylabel('Height')

outtab = np.vstack((zl,wavedrag[-1,:]/adimfac)).T
np.savetxt('wavedrag.dat',outtab)
