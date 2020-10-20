#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit

plt.ion()

#dir0 = '/home/bderembl/work/MITgcm/myrun/test_kw_energetics/run/'
dir0 = '/media/bderembl/workd/MITgcm/myrun/test_kw_energetics/run03/'
dir1 = dir0 + 'mnc*/'
dir2 = dir0 + 'mnc_test_0001/'

file0 = 'grid.t*'
file1 = 'state.*'
file2 = 'oceDiag.*'

alphat = 2e-4
go = 9.81

# grid
f0 = mit.mnc_files(dir1 + file0)

RC    = f0.variables['RC'][:]
DRC   = f0.variables['drC'][:]
DRF   = f0.variables['drF'][:]
RF    = f0.variables['RF'][:]

XC    = f0.variables['XC'][:,:]
YC    = f0.variables['YC'][:,:]

si_y,si_x = XC.shape
si_z      = RC.size

dx = XC[1,1] - XC[0,0]
dy = YC[1,1] - YC[0,0]
dz = RC[1] - RC[0]

dv = np.abs(dx*dy*dz)

f2 = mit.mnc_files(dir1 + file2)
T = f2.variables['T'][:]
si_t = len(T)

# compute KE, PE
circ = np.zeros((si_t))

for nt in range (0,si_t):
  rv = f2.variables['momVort3'][nt,0,:si_y,:si_x]

  circ[nt] = np.sum(np.sum(rv,0),0)


circ = circ/(si_y*si_x)

Td = T/86400

np.savetxt('mit_rv.dat',circ)


plt.figure()
plt.plot(T,circ[:],'k')

