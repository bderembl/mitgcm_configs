#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf

plt.ion()


dir0 = '../run/mnc_test_0004/'
file1 = 'state.0000000000.t001.nc'
file2 = 'grid.t001.nc'

f1 = netcdf.netcdf_file(dir0 + file1,'r')
f2 = netcdf.netcdf_file(dir0 + file2,'r')

topo = f2.variables['Depth'][:].copy().squeeze()

zl = f1.variables['Zl'][:].copy().squeeze()
xx = f1.variables['X'][:].copy().squeeze()

si_x = len(xx)
si_z = len(zl)

L0 = 1000
U0 = 1.0

zl = (zl - zl[-1])/L0
xx = (xx - (xx[0] + xx[-1])/2)/L0
topo = (topo[0] - topo)/L0

xg, zg = np.meshgrid(xx,zl)

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

plt.figure()
wmax = 0.1
wmin = -wmax
w = np.where(w>wmax,wmax,w)
w = np.where(w<wmin,wmin,w)
vcont = np.arange(-1,1,3e-3)
plt.contour(xx,zl,w,vcont,colors='k',linewidths=1)
plt.plot(xx,topo,linewidth=2)
plt.fill_between(xx, topo, 0*topo,facecolor="1",zorder=3)
plt.ylim([0, 1.0])

plt.figure()
vcont = np.arange(-1,1,1e-2)
plt.contour(xx,zl,ua,vcont,colors='k',linewidths=1)
plt.plot(xx,topo,linewidth=2)
plt.fill_between(xx, topo, 0*topo,facecolor="1",zorder=3)
plt.ylim([0, 1.0])

plt.figure()
nskx = 5
nskz = 1
nzm = np.argmin(np.abs(zl-0.2))
#nzm = np.argmin(np.abs(zl-0.6))

plt.plot(xx,topo,linewidth=2)
plt.quiver(xg[nzm::nskz,::nskx],zg[nzm::nskz,::nskx],u[nzm::nskz,::nskx],w[nzm::nskz,::nskx], angles='xy', scale_units='xy')
plt.fill_between(xx, topo, 0*topo,facecolor="1",zorder=3)
plt.xlim([-5, 5])

plt.figure()
vcont = np.linspace(0,10,80)
plt.contour(xx,zl,temp,vcont,colors='k',linewidths=1)
plt.plot(xx,topo,linewidth=2)
plt.fill_between(xx, topo, 0*topo,facecolor="1",zorder=3)
plt.xlim([-5, 5])
plt.ylim([0, 0.5])
