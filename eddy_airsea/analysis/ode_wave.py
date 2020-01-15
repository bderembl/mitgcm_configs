#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

plt.ion()

f0 = 1e-4
u0 = 1.0
R0 = 40e3      # radius
vmax = -1.0     # m/s



def v1(rr):
  v = -vmax*rr/R0*np.exp(-0.5*(rr/R0)**2)
#  v = -vmax*np.tanh(rr/R0)/(np.cosh(rr/R0))**2/(np.tanh(1.0)/(np.cosh(1.0))**2)
  return v
def dv1(rr):
  v = -vmax/R0*np.exp(-0.5*(rr/R0)**2)*(1-(rr/R0)**2)
#  v = -vmax*2/R0*np.tanh(rr/R0)/((np.cosh(rr/R0))**2)*(1/(np.cosh(rr/R0))**2 - (np.tanh(rr/R0))**2)/(np.tanh(1.0)/(np.cosh(1.0))**2)
  return v


def f(r, t):
  omega = np.sqrt((dv1(r)+v1(r)/r + f0)*(2*v1(r)/r + f0))
  return u0*np.sin(omega*t)

si_r = 30
si_t = 30000
r0 = np.linspace(1,5*R0,si_r)
t = np.linspace(0, si_t/f0/1000, si_t)
ra = np.zeros((si_t,si_r))
for ni in range(0,si_r):
  ra[:,ni] = integrate.odeint(f, r0[ni], t).squeeze()

plt.figure()
plt.plot(t*f0/(2*np.pi),ra/R0,'k',linewidth=1)
plt.xlabel(r'$tf/2\pi$')
plt.ylabel(r'$r_p/R_0$')

plt.xlim([np.min(t*f0/(2*np.pi)), np.max(t*f0/(2*np.pi))])
plt.ylim([np.min(ra/R0), 1.05*np.max(ra/R0)])

plt.savefig("ode_k0.pdf",bbox_inches='tight')
