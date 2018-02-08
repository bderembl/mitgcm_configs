#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

plt.ion()

flag_geom = 1 # cylindrical coordinates

# both rmax and r0 are non dimensional variables 
# the length unit is the deformation radius
rmax = 20  # tank radius
r0 = 10.0  #inner cylinder radius

si_r = 201 # number of points

rr = np.linspace(0,rmax,si_r)
dr = rr[1] - rr[0]

rr2 = 0.5*(rr[1:] + rr[:-1])


eta0 = np.zeros(si_r+1) 
eta0[1:-1] = np.where(rr2<r0,1,0) 
eta0[0] = eta0[1]
eta0[-1] = eta0[-2]
v0 = np.zeros(si_r) 

xx = np.concatenate((eta0,v0))


def f(x,si_r,rr,dr,eta0,flag_geom):
  eta = x[:si_r+1]
  v = x[si_r+1:]

  fout = np.zeros(2*si_r+1)

  detadr = (eta[1:] - eta[:-1])/(dr)
  detadr[0] = 0.
  detadr[-1] = 0.


  if flag_geom == 0:# cartesian
    fout[1:si_r] = (v[1:] - v[:-1])/(dr) - eta[1:-1] + eta0[1:-1]
    fout[si_r+1:] = -v + detadr
  if flag_geom == 1:# cylindrical
    fout[1:si_r] = (rr[1:]*v[1:] - rr[:-1]*v[:-1])/(dr*0.5*(rr[1:]+rr[:-1])) - eta[1:-1] + eta0[1:-1]
    fout[si_r+1:] = -v -v**2 + detadr

  return fout

sol = fsolve(f,xx,(si_r,rr,dr,eta0,flag_geom))
etaf = sol[1:si_r]
vf = sol[si_r+1:]

plt.figure()
plt.plot(rr2,etaf,'k',label='eta')
plt.plot(rr2,eta0[1:-1],'k--')
plt.plot(rr,vf,'r',label='v')
plt.xlabel('r/Rd')
plt.legend()
