#!/usr/bin/env python

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
import spoisson 
import def_radius
from scipy import interpolate
from scipy.interpolate import interp1d
import glob

plt.ion()

binprec = '>f4'

flag_conf = 2 # 0: samelson, 1: grooms 2: basilisk

#% ================== GRID =====================================
rSphere = 6370.e3
deg2m = 2*np.pi*rSphere/360.0
gg = 9.8
alphaT = 2e-4

si_x = 100
si_y = 100
#si_x = 720
#si_y = 720
if flag_conf == 0:
  si_z = 33
elif flag_conf == 1:
  si_z = 31
elif flag_conf == 2:
  si_z = 30


si_x1 = si_x + 1
si_y1 = si_y + 1

# in m
if flag_conf == 0 :
  Lx = 5000.0e3
  Ly = 5000.0e3
elif flag_conf == 1:
  Lx = 3000.0e3
  Ly = 3000.0e3
elif flag_conf == 2:
  Lx = 5000.0e3
  Ly = 5000.0e3

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

if flag_conf == 0:
  dz1 = np.array([  37.96964884,   64.27943707,   53.47713828,   55.25052547,
         57.14580417,   59.17549133,   61.35478616,   63.70082498,
         66.23372436,   68.97643209,   71.95606828,   75.20511746,
         78.76157761,   82.67134428,   86.99014783,   91.7853415 ,
         97.14066982,  103.16058993,  109.97712612,  117.75970459,
        126.72990561,  137.18292117,  149.52003956,  164.30348158,
        182.34416842,  204.85766232,  233.75503719,  272.22827077,
        326.05469227,  406.94121271,  543.09982806,  532.52164274,
        217.48963743])
elif flag_conf == 1:
  dz1 = np.zeros((si_z))
  dz1[0:4]   = 25.
  dz1[4:8]   = 50.
  dz1[8:12]  = 75.
  dz1[12:16] = 100.
  dz1[16:21] = 150.
  dz1[21:26] = 200.
  dz1[26:]   = 250.
elif flag_conf == 2:
  dz1 = 5000/si_z*np.ones((si_z))

# # 1 layer configuration
# si_z = 1
# dz1 = np.zeros((si_z))
# dz1[0] = 4000.0


zz = np.reshape(dz1.cumsum(),(si_z,1,1))

dz_fi2 = dz1/2.0
dz2 = np.array([dz_fi2[i] + dz_fi2[i+1] for i in range(len(dz_fi2)-1)])
dz2 = np.reshape(dz2[0:si_z-1],(si_z-1,1,1))


dx1.astype(binprec).tofile('dx.box')
dy1.astype(binprec).tofile('dy.box')
dz1.astype(binprec).tofile('dz.box')


# ==== physical parameters
if flag_conf == 0:
  fMin = 3.78e-05
  fMax = 1.3e-4
elif flag_conf == 1:
  fMin = 4.5e-5
  fMax = 1.0e-4
if flag_conf == 2:
  fMin = 3.e-05
  fMax = 1.3e-4

fmid = 0.5*(fMin + fMax)
  
beta = (fMax-fMin)/Ly
ff   = np.linspace(fMin,fMax,si_y)

print('f_south = {0}; beta = {1}'.format(fMin,beta) )


#%==================== LAND ===================================

landh  = np.zeros((si_y,si_x));


H = zz[-1].squeeze()
landh = -H + landh

# walls
landh[:,0] = 0.0
landh[-1,:] = 0.0
landh.astype(binprec).tofile('topog.box')

#%=============== Surface forcing ===================================

# -- temperature --
sst = np.zeros((si_y,si_x));

if flag_conf == 0:
  TS = 40.0 # should be 50, but I chose 40 # because the flux is a the top fo the ekman layer
  TN = 0.0
elif flag_conf == 1:
  TS = 22.0
  TN = 2.0
elif flag_conf == 2:
  TS = 22.0
  TN = 2.0

sst = (TN-TS)*yg/Ly + TS

#thetaClimFile
sst.astype(binprec).tofile('sstclim.box')

# relax time scale (grooms)
rho0 = 1023.0
Cp   = 4000.0
tauThetaClimRelax = rho0*Cp*dz1[0]/35. # 35 Watts per square meter per degree Kelvin

# relax time scale (samelson 97)
#(gamma*U*D/L/dz[1]) = 5*6e-6/37 ~ 15 days
#tauThetaClimRelax = 1233333.0
# I (arbitrarily..) set it to 50 days
tauThetaClimRelax = 4320000

# -- wind -- 
windx = np.zeros((si_y,si_x));

if flag_conf == 0:
  tauW = 0.4
elif flag_conf == 1:
  tauW = 0.2
elif flag_conf == 2:
  tauW = 0.4
windx = -tauW*np.sin(2*np.pi*yg/Ly )

windx = windx*ff.reshape(si_y,1)/fMin
windx.astype(binprec).tofile('windx.box')

#% ==============  background density profile ===================

# linear stratification
dep_l = np.linspace(0,H,si_z)

temp_f = (TN-TS)*(dep_l/H) + TS

if si_z > 1:
  # interpolate on the new vertical grid
  func2 = interp1d(dep_l, temp_f)
  temp_i = func2(zz)
else:
  temp_i = 1.0*temp_f

temp_i = temp_i.reshape((si_z,1,1))

temp_i.astype(binprec).tofile('tref.box')
#sref.astype(binprec).tofile('sref.box')


#%=============== initial conditions ===================================
# ### ideal ###

# uvel  = np.zeros((si_z,si_y,si_x));
# vvel  = np.zeros((si_z,si_y,si_x));
# theta = np.zeros((si_z,si_y,si_x));

# eta   = np.zeros((si_y,si_x));

# theta = theta + 4.0
# #theta = theta + temp_i - TN
# #theta = theta*(1-yg/Ly) + TN


# uvel.astype(binprec).tofile('uinit.box')
# vvel.astype(binprec).tofile('vinit.box')
# theta.astype(binprec).tofile('tinit.box')

# eta.astype(binprec).tofile('einit.box')


#### from PG ###

dir0 = './data_input/'
if flag_conf == 0:

  file1 = 'var_proj_s.nc'
  
  f1 = netcdf.netcdf_file(dir0 + file1,'r')
  
  uvel  = f1.variables['u' ][:,:,:].copy()
  vvel  = f1.variables['v' ][:,:,:].copy()
  theta = f1.variables['ti'][:,:,:].copy()
  
elif flag_conf == 2:
  # PG scales
  #L = 5000e3  # m
  H = 5000    # m
  beta = 2.0e-11 # 1/m/s
  N2 = 1e-6  #  (1/s**2)
  Bs = N2*H
  Thetas = Bs/gg/alphaT # 1/g alpha
  Us = N2*H**2/(beta*Lx**2)
  fnot = 3e-5
  gg = 9.80665 # nemo value
  
  ff = fnot + beta*yg # should be at u and v points
  fmid = fnot + 0.5*Ly*beta  
  
  fileb = 'b*'
  fileu = 'u*'
  
  allfilesb = sorted(glob.glob(dir0 + fileb));
  allfilesu = sorted(glob.glob(dir0 + fileu));
  
  # dimensions
  b = np.fromfile(allfilesb[0],'f4')
  N = int(b[0])
  N1 = N + 1
  nl2 = int(len(b)/N1**2)
  nl = nl2 - 2
  
  b = np.fromfile(allfilesb[-1],'f4').reshape(nl2,N1,N1).transpose(0,2,1)
  uv = np.fromfile(allfilesu[-1],'f4').reshape(2*nl2,N1,N1).transpose(0,2,1)
  
  theta = Thetas*(b[1:-1,1:,1:] - b.min()) + 2.0
  uvel = Us*uv[2:-2:2,1:,1:]
  vvel = Us*uv[3:-2:2,1:,1:]

  
si_zpg,si_ypg,si_xpg = theta.shape
dxpg = dx*si_x/si_xpg

# compute pressure for SSH
dudy = np.diff(uvel,1,1)/dxpg
dvdx = np.diff(vvel,1,2)/dxpg

vort = dvdx[0,:-1,:] - dudy[0,:,:-1]

psi = spoisson.sol(vort[:])
psi = psi.reshape((si_xpg-1,si_ypg-1))
psi = psi*dxpg*dxpg*fmid/gg

eta = np.zeros((si_ypg,si_xpg))
eta[:-1,:-1] = psi

# old grid
xx = np.linspace(0,1,si_xpg)
yy = np.linspace(0,1,si_ypg)

xog,yog = np.meshgrid(xx,yy)



xn = np.linspace(0,1,si_x)
yn = np.linspace(0,1,si_y)

xng,yng = np.meshgrid(xn,yn)


uvel_n  = np.zeros((si_z,si_y,si_x))
vvel_n  = np.zeros((si_z,si_y,si_x))
theta_n = np.zeros((si_z,si_y,si_x))

eta_n = np.zeros((si_y,si_x))

for nz in range(0,si_z):
  fint = interpolate.interp2d(xx, yy,uvel[nz,:,:], kind='cubic')
  uvel_n[nz,:,:] = fint(xn,yn)
  
  fint = interpolate.interp2d(xx, yy,vvel[nz,:,:], kind='cubic')
  vvel_n[nz,:,:] = fint(xn,yn)

  fint = interpolate.interp2d(xx, yy,theta[nz,:,:], kind='cubic')
  theta_n[nz,:,:] = fint(xn,yn)

fint = interpolate.interp2d(xx, yy,eta, kind='cubic')
eta_n = fint(xn,yn)

#np.savetxt('upg.dat',uvel_n[0,:,:])
#np.savetxt('sstpg.dat',theta_n[0,:,:])

uvel_n.astype(binprec).tofile('uinit.box')
vvel_n.astype(binprec).tofile('vinit.box')
theta_n.astype(binprec).tofile('tinit.box')

eta_n.astype(binprec).tofile('einit.box')

#---------------------
# ------ RBCS --------
#---------------------

tmask  = np.ones((si_z,si_y,si_x))
tmask.astype(binprec).tofile('tmask.box')

# compute relaxation length scale

N2 = -gg*alphaT*np.diff(theta_n,axis=0)/dz2
N2_min = 1e-7
N2 = np.where(N2<N2_min, N2_min, N2)

gp = N2*dz2
lmax = Lx/10
filt_len = np.zeros((si_y,si_x))
for nx in range(0,si_x):
  for ny in range(0,si_y):
    rd = def_radius.cal_rad(dz1,gp[:,ny,nx],ff[ny,nx])
    filt_len[ny,nx] = np.min([10*rd[1],lmax])

filt_len.astype(binprec).tofile('filter_length.box')
