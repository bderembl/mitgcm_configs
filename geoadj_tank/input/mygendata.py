#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt

plt.ion()

binprec = '>f4'


H = 0.15      # water level in the tank
H1 = 0.05      # thickness of the lens
R0 = 0.08      # inner cylinder radius
R1 = 0.5       # radius of the tank
Omega = 1.7453 # rotation rate
S0 = 5.0       # upper layer salinity
S1 = 30.0      # lower layer salinity

f0 = 2*Omega


# ========== grid =========

si_x = 100
si_y = 100
si_z = 29

Lx = 2*R1
Ly = 2*R1

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

x_c = Lx/2
y_c = Ly/2

rad_gg = np.sqrt((xg-x_c)**2 + (yg-y_c)**2)
rad_cc = np.sqrt((xc-x_c)**2 + (yc-y_c)**2)
rad_gu = np.sqrt((xu-x_c)**2 + (yu-y_c)**2)
rad_gv = np.sqrt((xv-x_c)**2 + (yv-y_c)**2)

theta_gg = np.arctan2(yg-y_c,xg-x_c)
theta_cc = np.arctan2(yc-y_c,xc-x_c)
theta_gu = np.arctan2(yu-y_c,xu-x_c)
theta_gv = np.arctan2(yv-y_c,xv-x_c)


# ============= topography =============

h_mit = np.where(rad_gg>R1,0.0,-H)
h_mit[0,:] = 0
h_mit[:,0] = 0

h_mit.astype(binprec).tofile('topo.box')


# ========= initial conditions ==========

sinit = S1 + np.zeros((si_z,si_y,si_x))

for k in range(0,si_z):
  if (zz[k]<H1):
    sinit[k,:,:] = np.where(rad_gg<R0,S0,S1)

uinit = 0.01*(0.5-np.random.rand(si_z,si_y,si_x))
vinit = 0.01*(0.5-np.random.rand(si_z,si_y,si_x))

uinit.astype(binprec).tofile('uinit.box')
vinit.astype(binprec).tofile('vinit.box')
sinit.astype(binprec).tofile('sinit.box')

