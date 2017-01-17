#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.special, scipy.interpolate, scipy.integrate
import scipy.optimize
from scipy import special

plt.ion()

binprec = '>f4'

#% ================ CONST =========================================                                  

H = -500       # mean depth
fo = 2.53e-5   # coriolis parameter NBC 4*np.pi/86400*np.sin(10*np.pi/180)
beta = 2.2e-11 # 4*np.pi/86400*np.cos(10*np.pi/180)/6400e3
g = 9.81       # gravity   

# size of the domain
Lx = 500.0e3   # (in m)
Ly = 500.0e3   # (in m)


#% ================ GRID =========================================                                  

# number of grid points
si_x = 400; 
si_y = 400;
si_z = 1;

si_x1 = si_x + 1
si_y1 = si_y + 1

nx1 = 1;
nx2 = np.floor(si_x/2);
ny1 = 1;
ny2 = np.floor(si_y/2);

dx = Lx/si_x;
dy = Ly/si_y;
dz = -H



xx = Lx*(np.arange(0,si_x) + 0.5)/(1.0*si_x)
yy = Ly*(np.arange(0,si_y) + 0.5)/(1.0*si_y)

xx1 = Lx*(np.arange(0,si_x+1) )/(1.0*si_x)
yy1 = Ly*(np.arange(0,si_y+1) )/(1.0*si_y)


xg,yg = np.meshgrid(xx,yy) 
xu,yu = np.meshgrid(xx1[:-1],yy) 
xv,yv = np.meshgrid(xx,yy1[:-1]) 
xc,yc = np.meshgrid(xx1,yy1) 


dx1 = dx*np.ones((si_x))
dy1 = dy*np.ones((si_y))
dz1 = dz*np.ones((si_z))

dx1.astype(binprec).tofile('dx.box')
dy1.astype(binprec).tofile('dy.box')
dz1.astype(binprec).tofile('dz.box')

#% ================ Topography =============================  


h_mit = H*np.ones((si_y,si_x))

h_mit[0,:] = 0.0
h_mit[:,0] = 0.0

# corner:
h_mit[0:int(si_y/2),0:int(si_x*2/5)] = 0.0

h_mit.astype(binprec).tofile('topo.box')


#% ================ initial condition =============================  

# Initial position of the vortex:
x_c = 250e3     # center (x, km)
y_c = 150e3     # center (y, km)
R0 = 20e3       # radius (km)

Gamma = -2.0e4   # strength of the vortex

rad_gg = np.sqrt((xg-x_c)**2 + (yg-y_c)**2)
rad_cc = np.sqrt((xc-x_c)**2 + (yc-y_c)**2)
rad_gu = np.sqrt((xu-x_c)**2 + (yu-y_c)**2)
rad_gv = np.sqrt((xv-x_c)**2 + (yv-y_c)**2)

theta_gg = np.arctan2(yg-y_c,xg-x_c)
theta_cc = np.arctan2(yc-y_c,xc-x_c)
theta_gu = np.arctan2(yu-y_c,xu-x_c)
theta_gv = np.arctan2(yv-y_c,xv-x_c)

u1 = Gamma*rad_gu/(2*np.pi*R0**2)*np.sin(-theta_gu)
u2 = Gamma/(2*np.pi*rad_gu)*np.sin(-theta_gu)

v1 = Gamma*rad_gv/(2*np.pi*R0**2)*np.cos(theta_gv)
v2 = Gamma/(2*np.pi*rad_gv)*np.cos(theta_gv)

u_out = np.where(rad_gu>R0,u2,u1)
v_out = np.where(rad_gv>R0,v2,v1)


# lamb-oseen vortex
u_out = Gamma/(2*np.pi*rad_gu)*(1-np.exp(-rad_gu**2/R0**2))*np.sin(-theta_gu)
v_out = Gamma/(2*np.pi*rad_gv)*(1-np.exp(-rad_gv**2/R0**2))*np.cos(theta_gv)

vort =  Gamma/(np.pi*R0**2)*(np.exp(-rad_cc**2/R0**2))
psi = 0.0*vort

dudx = np.diff(u_out)/dx
dvdy = np.diff(v_out,1,0)/dy

div = dudx[:si_y-1,:si_x-1] + dvdy[:si_y-1,:si_x-1]



rr = np.linspace(0.0,2*Lx, 10*si_x)

# orthoradial velocity
# vtheta = Gamma/(2*np.pi*rr)*(1-np.exp(-rr**2/R0**2))

def integrand(rr):
  # cyclostrophic balance (f = 0)
  res = (Gamma/(2*np.pi*rr)*(1-np.exp(-rr**2/R0**2)))**2/rr 
  # geostrophic balance
  res = res + fo*(Gamma/(2*np.pi*rr)*(1-np.exp(-rr**2/R0**2)))
  res = np.where(rr == 0, 0.0,res)
  return res

def comp_p(x):
  
  if x ==0:
    xx = 1e-12
  else:
    xx = 1.0*x

  a,b = scipy.integrate.quad(integrand,0,xx)
  return a

# compute pressure
p1 = [ comp_p(x) for x in rr.flatten() ]
fint = scipy.interpolate.interp1d(rr, p1)

# derive sea surface elevation
e_out = 1/(g)*(Gamma/(2*np.pi))**2* 0.5*(scipy.special.expi(-rad_gg**2/R0**2)/R0**2 - (1-np.exp(-rad_gv**2/R0**2))/rad_gg**2)

e_out2 = fint(rad_gg)/g
u_2 = Gamma/(2*np.pi*rad_gg)*(1-np.exp(-rad_gg**2/R0**2))*np.cos(theta_gg)
v_2 = Gamma/(2*np.pi*rad_gg)*(1-np.exp(-rad_gg**2/R0**2))*np.sin(-theta_gg)

dudx2,dudy2 = np.gradient(u_2)
dvdx2,dvdy2 = np.gradient(v_2)

vort2 = dvdx2 - dudy2
div2 = dudx2 + dvdy2



dvdx = np.diff(v_out,1,1)/dx
dudy = np.diff(u_out,1,0)/dy

vort_diag = dvdx[:si_y-1,:] - dudy[:,:si_x-1]

u_out.astype(binprec).tofile('uinit.box')
v_out.astype(binprec).tofile('vinit.box')

e_out2.astype(binprec).tofile('einit.box')
