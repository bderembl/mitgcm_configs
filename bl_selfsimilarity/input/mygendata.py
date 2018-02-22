#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt

plt.ion()

binprec = '>f4'


Lx = 1.0
Ly = 1.0
H = 100.0
U0 = 1.0

# ========== grid =========

si_x = 1
si_y = 1
si_z = 1000

dx = Lx/si_x
dy = Ly/si_y
dz = H/si_z

dx1 = dx*np.ones((si_x))
dy1 = dy*np.ones((si_y))
dz1 = dz*np.ones((si_z))

xx = Lx*(np.arange(0,si_x) + 0.5)/(1.0*si_x)
yy = Ly*(np.arange(0,si_y) + 0.5)/(1.0*si_y)

xx1 = Lx*(np.arange(0,si_x+1) )/(1.0*si_x)
yy1 = Ly*(np.arange(0,si_y+1) )/(1.0*si_y)

zz = np.cumsum(dz1)

xg,yg = np.meshgrid(xx,yy) 
xu,yu = np.meshgrid(xx1[:-1],yy) 
xv,yv = np.meshgrid(xx,yy1[:-1]) 
xc,yc = np.meshgrid(xx1,yy1) 

dx1.astype(binprec).tofile('dx.box')
dy1.astype(binprec).tofile('dy.box')
dz1.astype(binprec).tofile('dz.box')


# ============= topography =============

h_mit = -H + np.zeros((si_y,si_x))

h_mit.astype(binprec).tofile('topo.box')



# ========= initial conditions ==========

uinit = U0 + np.zeros((si_z,si_y,si_x))

uinit.astype(binprec).tofile('uinit.box')

