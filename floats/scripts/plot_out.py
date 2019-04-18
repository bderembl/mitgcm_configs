#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit

plt.ion()

dir1 = '../run/'

file_u = 'U*'
file_v = 'V*'
file_t = 'T*'
file_s = 'S*'
file_e = 'Eta*'

iters = mit.mds.scanforfiles(dir1 + file_u)
si_t = len(iters)

# grid files
XC = mit.rdmds(dir1 + 'XC').squeeze()
YC = mit.rdmds(dir1 + 'YC').squeeze()


nt = -1

U = mit.rdmds(dir1 + file_u, itrs = [iters[nt]]).squeeze()
V = mit.rdmds(dir1 + file_v, itrs = [iters[nt]]).squeeze()
T = mit.rdmds(dir1 + file_t, itrs = [iters[nt]]).squeeze()
S = mit.rdmds(dir1 + file_s, itrs = [iters[nt]]).squeeze()
e = mit.rdmds(dir1 + file_e, itrs = [iters[nt]]).squeeze()


plt.figure()
plt.contourf(XC*1e-3,YC*1e-3,S[0,:,:],20)
#plt.contourf(XC*1e-3,YC*1e-3,V[0,:,:],20)
plt.colorbar()
plt.contour(XC*1e-3,YC*1e-3,e,colors='k')

plt.xlabel('x (km)')
plt.ylabel('y (km)')
