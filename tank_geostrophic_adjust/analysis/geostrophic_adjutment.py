#!/usr/bin/env python

"""This program finds the steady state solution of the geostrophic
adjustment the equations are: (i) the geostrophic balance, and (ii)
the conservation of potential vorticity. In cartesian geometry, these
two equations are

(i)    -v = dh/dx
(ii)   dv/dx - h = -h_0

with v the velocity and h the elevation of the interface

the program returns the solution vector (v,h) and plots the initial
condition (h_0 dashed) and final state (h)

v and h are discretized on a C-grid where the velocity is at the face
of the cell and the height of the interface is at the center of the
cell.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

plt.ion()

flag_geom = 1 # geometry: 0: cartesian, 1: cylindrical coordinates

# both rmax and r0 are non dimensional variables 
# the length unit is the deformation radius
rmax = 20  # tank radius
r0 = 10.0   #inner cylinder radius

Nx = 100     # number of points
si_r = 2*Nx + 1 

rr = np.linspace(0,rmax,si_r)
dr = rr[1] - rr[0]

rr2 = 0.5*(rr[1:] + rr[:-1])

# initial condition
h0 = np.zeros(si_r+1) 
h0[1:-1] = np.where(rr2<r0,-1,0) 
h0[0] = h0[1]
h0[-1] = h0[-2]
v0 = np.zeros(si_r) 

xx = np.concatenate((h0,v0))


def f(x,si_r,rr,dr,h0,flag_geom):
  h = x[:si_r+1]
  v = x[si_r+1:]

  fout = np.zeros(2*si_r+1)

  dhdr = (h[1:] - h[:-1])/dr
  dhdr[0] = 0.
  dhdr[-1] = 0.


  if flag_geom == 0:# cartesian
    fout[1:si_r] = (v[1:] - v[:-1])/(dr) - h[1:-1] + h0[1:-1]
    fout[si_r+1:] = -v + dhdr
  if flag_geom == 1:# cylindrical
    fout[1:si_r] = (rr[1:]*v[1:] - rr[:-1]*v[:-1])/(dr*0.5*(rr[1:]+rr[:-1])) - h[1:-1] + h0[1:-1]
    fout[si_r+1:] = -v -v**2/rr + dhdr
    fout[si_r+1] = 0.0

  return fout

sol = fsolve(f,xx,(si_r,rr,dr,h0,flag_geom))
hf = sol[1:si_r]
vf = sol[si_r+1:]

plt.figure()
plt.plot(rr2,hf,'k',label='h',linewidth=1)
plt.plot(rr2,h0[1:-1],'k--',linewidth=1)
#plt.plot(rr,vf,'r',label='v',linewidth=1)
plt.xlabel('r/Rd')
plt.ylabel('h/h0')
plt.legend()
