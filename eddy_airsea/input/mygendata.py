#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate, scipy.integrate

plt.ion()

binprec = '>f4'

flag_plot = 0
flag_bt = 0  # 1: barotropic vortex, 0: baroclinic vortex

# # physical constants
rho_const = 999.8
alphaK = 2.0e-4
g0 = 9.8
f0 = 1e-4


#% ================== NEW GRID =====================================

def stretch(xf,yf,Lx,si_x,rev):

  hh = np.linspace(0,1,si_x+1)
  xf = Lx*np.interp(hh,xf,yf)

  dx = np.diff(xf)

  # reverse order to get high resolution near the bottom
  if rev:
    dx = dx[::-1]
  xf[1:] = np.cumsum(dx)
  
  xc = xf[0:-1] + 0.5*dx
  return xc,xf,dx

si_x = 250
si_y = 250
si_z = 20

si_x1 = si_x + 1
si_y1 = si_y + 1
si_z1 = si_z + 1

# in m
Lx = 300.0e3
Ly = 300.0e3
Lz = 1300

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

slope = 4
xfz = np.linspace(0,1,1000)
yfz = np.sinh(slope*xfz)/np.sinh(slope)
zc,zf,dz1 = stretch(xfz,yfz,Lz,si_z,0)


iz = np.argmin(np.abs(zf-500.0))

print ('dx= ', dx)
print ('min dz: ', np.min(dz1))
print ('max dz: ', np.max(dz1))
print ('nb layers above 500m:', iz, '/', si_z)

if np.sum(dz1 < 0) > 0:
  print ('you need you change the polynomial fit!')


dx1.astype(binprec).tofile('dx.box')
dy1.astype(binprec).tofile('dy.box')
dz1.astype(binprec).tofile('dz.box')


#% ============== background density profile ===================
N2 = 3e-5
#N2 = 0.
tref = -N2/g0/alphaK*zc
tref = tref - tref[-1]

# tref_dat = np.loadtxt('temp_winter.dat')
# # add upper and lower boundary for interpolation
# si_tref,naux = tref_dat.shape
# tref_dat2 = np.zeros((si_tref+2,2))
# tref_dat2[1:-1,:] = tref_dat
# tref_dat2[0,1] = tref_dat[0,1]
# tref_dat2[-1,:] = tref_dat[-1,:]
# tref_dat2[-1,0] = -10000

# t_interp = scipy.interpolate.interp1d(tref_dat2[:,0], tref_dat2[:,1])
# tref = t_interp(-zc)

tref = tref.reshape((si_z,1,1))

#%==================== SST - LAND ===================================

landh  = np.zeros((si_y,si_x));
theta  = np.zeros((si_z,si_y,si_x));
uvel   = np.zeros((si_z,si_y,si_x));
vvel   = np.zeros((si_z,si_y,si_x));


H = dz1.cumsum()[-1]
landh = -H + landh

#gaussian eddy
x_c = Lx/2     # x center
y_c = Ly/2     # y center
R0 = 40e3      # radius

z0 = 200.0     # characteristic depth (m)
vmax = 0.5     # m/s

# vertical profile
FZ = 1-scipy.special.erf(zc/z0)
if flag_bt:
  FZ = 0*FZ + 1
FZ = FZ.reshape(si_z,1,1)

# vertical derivative
FpZ = -1/z0*2/np.sqrt(np.pi)*np.exp(-zc**2/z0**2)
if flag_bt:
  FpZ = 0*FpZ 
FpZ = FpZ.reshape(si_z,1,1)

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
def geostrophic_part(rr):
  v = vel_rankine(rr)
  res = f0*v
  res = np.where(rr == 0, 0.0,res)
  return res

def cyclostrophic_part(rr):
  v = vel_rankine(rr)
  res = v**2/rr
  res = np.where(rr == 0, 0.0,res)
  return res

def comp_p(x, func):
  if x ==0:
    xx = 1e-12
  else:
    xx = 1.0*x

  a,b = scipy.integrate.quad(func,0,xx)
  return a

rr = np.linspace(0.0,2*Lx, 10*si_x)
p1 = [ comp_p(x,geostrophic_part) for x in rr.flatten() ]
p2 = [ comp_p(x,cyclostrophic_part) for x in rr.flatten() ]
fint1 = scipy.interpolate.interp1d(rr, p1)
fint2 = scipy.interpolate.interp1d(rr, p2)

p_out1 = rho_const*fint1(rad_gg)
p_out2 = rho_const*fint2(rad_gg)
# remove const at infinity
p_out1 = p_out1 - p_out1[0,0]
p_out2 = p_out2 - p_out2[0,0]

dpdz = FpZ*np.tile(p_out1,[si_z,1,1]) +  2*FZ*FpZ*np.tile(p_out2,[si_z,1,1])
rhop = dpdz/g0

# convert to temperature
theta_a = -rhop/(rho_const*alphaK) 
theta = theta_a + tref

# free surface
pres = FZ*np.tile(p_out1,[si_z,1,1]) + FZ**2*np.tile(p_out2,[si_z,1,1])
eta = pres[0,:,:]/rho_const/g0 - rhop[0,:,:]/rho_const*0.5*dz1[0]

uvel[0,:,:] = 0*uvel[0,:,:]
vvel[0,:,:] = 0*vvel[0,:,:]

uvel.astype(binprec).tofile('uinit.box')
vvel.astype(binprec).tofile('vinit.box')
theta.astype(binprec).tofile('tinit.box')

eta.astype(binprec).tofile('einit.box')

tref.astype(binprec).tofile('tref.box')
#sref.astype(binprec).tofile('sref.box')


landh.astype(binprec).tofile('topog.box')

# # === RBCS files
tmask  = np.zeros((si_z,si_y,si_x));
tmask[-1,:,:] = 1.0
trelax = 1.0*theta

tmask.astype(binprec).tofile('tmask.box')
trelax.astype(binprec).tofile('trelax.box')

# ==== atmospheric files
u0 = 10  # m/s
v0 = 10  # m/s
s0 = 240
t2 = 15
q2 = 1e-3

uatm  = u0*np.ones((si_y,si_x));
vatm  = v0*np.ones((si_y,si_x));
solar = s0*np.ones((si_y,si_x));
tair  = t2*np.ones((si_y,si_x));
qair  = q2*np.ones((si_y,si_x));

uatm.astype(binprec).tofile('u10.box')
vatm.astype(binprec).tofile('v10.box')
solar.astype(binprec).tofile('ssrd.box')
tair.astype(binprec).tofile('t2.box')
qair.astype(binprec).tofile('d2.box')


#====== figures ========
vmin = np.min(vvel)
vcont = np.linspace(vmin,-vmin,6)

plt.figure()
plt.contourf(xx*1e-3,-zc,theta[:,int(si_y/2),:],20)
plt.colorbar(label="Theta (K)")
plt.contour(xx*1e-3,-zc,vvel[:,int(si_y/2),:],vcont,colors='k')
plt.xlabel("x (km)")
plt.ylabel("z (m)")
plt.title("Temperature (color) and velocity (contours)")
plt.savefig("surf_eddy.png")
