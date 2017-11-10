#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.special, scipy.interpolate, scipy.integrate

plt.ion()

binprec = '>f8'

#% ================== NEW GRID =====================================

si_x = 100
si_y = 100
si_z = 40


si_x1 = si_x + 1
si_y1 = si_y + 1

# in m
Lx = 100.0e3
Ly = 100.0e3

dx = Lx/si_x;
dy = Ly/si_y;

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


dz_fi = 1.0*np.array([5,5,5,5,10,10,10,10,15,15,15,15,20,20,20,20,30,30,30,30,40,40,40,50,50,50,60,60,60,80,80,100,100,120,120,140,140,160,160,180,180,180,180,200,200,200,200,200,200,200,200])

dz_fi = dz_fi[::-1]

dz_fi2 = dz_fi/2.0

dz1 = dz_fi[0:si_z]

zz = np.reshape(dz1.cumsum(),(si_z,1,1))

dz2 = np.array([dz_fi2[i] + dz_fi2[i+1] for i in range(len(dz_fi2)-1)])
dz2 = np.reshape(dz2[0:si_z-1],(si_z-1,1,1))


dx1.astype(binprec).tofile('dx.box')
dy1.astype(binprec).tofile('dy.box')
dz1.astype(binprec).tofile('dz.box')

#%==================== dynamic variables ===================================

topog  = np.zeros((si_y,si_x));
theta  = np.zeros((si_z,si_y,si_x));
uvel   = np.zeros((si_z,si_y,si_x));
vvel   = np.zeros((si_z,si_y,si_x));


H = dz1.cumsum()[-1]
# flat topo
topog = -H + topog

lt = 20e3 # topo wave length
Ht = 100 # topo max height
topog = topog + Ht*(1+np.sin(2*np.pi*xg/lt)*np.sin(2*np.pi*yg/lt))

# # physical constants
rho_const = 999.8
alphaK = 2.0e-4
g0 = 9.8
f0 = 1e-4

# background density profile
N2 = 1e-6
temp_i = -N2/g0/alphaK*zz
temp_i = temp_i - temp_i[-1]

# eddy parameters
x_c = Lx/2
y_c = Ly/2
DeltaT = 5.0
t0 = temp_i[0] 
R0 = 14e3
vmax = 0.1

# vertical profile

# # surface intensified
# z0 = 600.0 # m
# FZ = np.exp(-zz**2/z0**2)
# FZ = FZ.reshape(si_z,1,1)

# # vertical derivative
# FpZ = -2*zz/z0**2*np.exp(-zz**2/z0**2)
# FpZ = FpZ.reshape(si_z,1,1)


# barotropic
FZ  = 0*zz + 1.0
FpZ = 0.0*zz


# grid at U,V,T points
rad_gg = np.sqrt((xg-x_c)**2 + (yg-y_c)**2)
rad_cc = np.sqrt((xc-x_c)**2 + (yc-y_c)**2)
rad_gu = np.sqrt((xu-x_c)**2 + (yu-y_c)**2)
rad_gv = np.sqrt((xv-x_c)**2 + (yv-y_c)**2)

theta_gg = np.arctan2(yg-y_c,xg-x_c)
theta_cc = np.arctan2(yc-y_c,xc-x_c)
theta_gu = np.arctan2(yu-y_c,xu-x_c)
theta_gv = np.arctan2(yv-y_c,xv-x_c)

# create velocity profile
# hyperpolic vortex (~ rankine)
def vel_rankine(rr):
  v = -vmax*np.tanh(rr/R0)/(np.cosh(rr/R0))**2/(np.tanh(1.0)/(np.cosh(1.0))**2)
  v = np.where(rr == 0, 0.0,v)
  return v

u_out = vel_rankine(rad_gu)*np.sin(-theta_gu)
v_out = vel_rankine(rad_gv)*np.cos(theta_gv)


# 3D velocity field
uvel = FZ*np.tile(u_out,[si_z,1,1])
vvel = FZ*np.tile(v_out,[si_z,1,1])


# compute pressure field
def integrand(rr):
  v = vel_rankine(rr)
  res = v**2/rr  + f0*v
  # no cyclostrophic
#  res = f0*v
  res = np.where(rr == 0, 0.0,res)
  return res

def comp_p(x):
  
  if x ==0:
    xx = 1e-12
  else:
    xx = 1.0*x

  a,b = scipy.integrate.quad(integrand,0,xx)
  return a

rr = np.linspace(0.0,2*Lx, 10*si_x)
p1 = [ comp_p(x) for x in rr.flatten() ]
fint = scipy.interpolate.interp1d(rr, p1)

p_out = rho_const*fint(rad_gg)
# remove const at infinity
p_out = p_out - p_out[0,0]


eta = p_out/rho_const/g0
pres = FZ*np.tile(p_out,[si_z,1,1])



# minus sign
#dpdz = np.diff(pres,1,0)/dz2
dpdz = FpZ*np.tile(p_out,[si_z,1,1])
rhop = dpdz/g0

# convert to temperature
theta_a = -rhop/(rho_const*alphaK) 
theta = theta_a + temp_i

uvel.astype(binprec).tofile('uinit.box')
vvel.astype(binprec).tofile('vinit.box')
theta.astype(binprec).tofile('tinit.box')

eta.astype(binprec).tofile('einit.box')

temp_i.astype(binprec).tofile('tref.box')

topog.astype(binprec).tofile('topog.box')

########## OBCS files ================

# East
u_e = uvel[:,:,-1]
v_e = vvel[:,:,-1]
t_e = theta[:,:,-1]

u_e.astype(binprec).tofile('u_E.box')
v_e.astype(binprec).tofile('v_E.box')
t_e.astype(binprec).tofile('t_E.box')

# East
u_w = uvel[:,:,0]
v_w = vvel[:,:,0]
t_w = theta[:,:,0]

u_w.astype(binprec).tofile('u_W.box')
v_w.astype(binprec).tofile('v_W.box')
t_w.astype(binprec).tofile('t_W.box')

# North
u_n = uvel[:,-1,:]
v_n = vvel[:,-1,:]
t_n = theta[:,-1,:]

u_n.astype(binprec).tofile('u_N.box')
v_n.astype(binprec).tofile('v_N.box')
t_n.astype(binprec).tofile('t_N.box')

# South
u_s = uvel[:,0,:]
v_s = vvel[:,0,:]
t_s = theta[:,0,:]

u_s.astype(binprec).tofile('u_S.box')
v_s.astype(binprec).tofile('v_S.box')
t_s.astype(binprec).tofile('t_S.box')

