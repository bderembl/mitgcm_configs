#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate, scipy.integrate

import re

plt.ion()

binprec = '>f4'


Lz = 0.1        # water level in the tank
R1 = 0.5       # radius of the tank
Omega = 1.7453 # rotation rate (2pi/T)
delta_omega = 0.1
# vertical grid stretching parameter
stretch_coef = 4

rho_const = 999.8
g0 = 9.8

f0 = 2*Omega

def stretch(xf,yf,Lx,si_x,rev=0):

  hh = np.linspace(0,1,si_x+1)
  xf = Lx*np.interp(hh,xf,yf)

  dx = np.diff(xf)

  # reverse order to get high resolution near the bottom
  if rev:
    dx = dx[::-1]
  xf[1:] = np.cumsum(dx)
  
  xc = xf[0:-1] + 0.5*dx
  return xc,xf,dx

# ========== grid =========

si_x = 64
si_y = 64
si_z = 29

Lx = 2*R1
Ly = 2*R1

dx = Lx/si_x
dy = Ly/si_y
dz = Lz/si_z

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


# vertical grid
xfz = np.linspace(0,1,1000)
yfz = np.sinh(stretch_coef*xfz)/np.sinh(stretch_coef)

zc,zf,dz1 = stretch(xfz,yfz,Lz,si_z,rev=1)

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

# number of processes
npx = int(max(si_x/64,1))
npy = int(max(si_x/64,1))


# ============= topography =============

h_mit = np.where(rad_gg>R1,0.0,-Lz)
h_mit[0,:] = 0
h_mit[:,0] = 0

h_mit.astype(binprec).tofile('topo.box')


# ========= initial conditions ==========

def vel_profile(rr):
  v = delta_omega*rr
  v = np.where(rr == 0, 0.0,v)
  return v

u_out = vel_profile(rad_gu)*np.sin(-theta_gu)
v_out = vel_profile(rad_gv)*np.cos(theta_gv)

# 3D velocity field
uvel = np.tile(u_out,[si_z,1,1])
vvel = np.tile(v_out,[si_z,1,1])


# compute pressure field
def integrand(rr):
  v = vel_profile(rr)

  res = v**2/rr  + f0*v
  res = np.where(rr == 0, 0.0,res)
  return res

def comp_p(x):
  
  if x ==0:
    xx = 1e-12
  else:
    xx = 1.0*x

  a,b = scipy.integrate.quad(integrand,0,xx)
  return a

Lmax = np.max([Lx,Ly])
si_max = np.max([si_x,si_y])
rr = np.linspace(0.0,2*Lmax, 10*si_max)
p1 = [ comp_p(x) for x in rr.flatten() ]
fint = scipy.interpolate.interp1d(rr, p1)

p_out   =  rho_const*fint(rad_gg)

# remove const at infinity
p_out = p_out - p_out[0,0]

eta = p_out/rho_const/g0


uvel.astype(binprec).tofile('uinit.box')
vvel.astype(binprec).tofile('vinit.box')
eta.astype(binprec).tofile('einit.box')

#=================== Write data files  =====================

with open("../code/SIZE.h_0", "r") as sources:
    lines = sources.readlines()
with open("../code/SIZE.h", "w") as sources:
    for line in lines:
      line2 = re.sub(r'MYSNX', str(int(si_x/npx)), line)
      line2 = re.sub(r'MYSNY', str(int(si_y/npy)), line2)
      line2 = re.sub(r'MYNPX', str(npx), line2)
      line2 = re.sub(r'MYNPY', str(npy), line2)
      line2 = re.sub(r'MYNR', str(si_z), line2)
      sources.write(line2)

