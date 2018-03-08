#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf

plt.ion()

# directory of the experiment
dir0 = '../run/mnc_test_0001/'
file1 = 'state.0000000000.t001.nc'

f1 = netcdf.netcdf_file(dir0 + file1, 'r')

# load coordinate: X, Z in (m), T in (s)
X = f1.variables['X'][:].copy()
Z = f1.variables['Z'][:].copy()
T = f1.variables['T'][:].copy()

si_x = len(X)
si_z = len(Z)
si_t = len(T)

# load variables: S (psu), U, W (m/s)
S = f1.variables['S'][:,:,:,:].copy().squeeze()
U = f1.variables['U'][:,:,:,:-1].copy().squeeze()
W = f1.variables['W'][:,:,:,:].copy()

# select time frame index for snapshot. actual time is T[nt]
nt = 10

plt.figure()
plt.contourf(X,Z,S[nt,:,:])
plt.axis('scaled')
plt.xlabel('X (m)')
plt.ylabel('Z (m)')

# find front position: -1 is for lowermost grid point
# we skip the first grid point because it is the wall
dS = np.diff(S[:,-1,1:],1,1)

# find index of the front
Sf_i = np.argmax(dS,axis=1)

# we plot x_f - x[end] (since the front propagates from right to left)
plt.figure()
plt.plot(T, X[-1] - X[Sf_i])
plt.xlabel('Time (s)')
plt.ylabel('Distance from init (m)')

