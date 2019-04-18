
#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

plt.ion()

binprec = '>f4'
flag_plot = 0


#% ================== GRID =====================================
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



rSphere = 6370.e3
deg2m = 2*np.pi*rSphere/360.0
gg = 9.8

Lx = 5000.0e3
Ly = 5000.0e3
Lz = 5000.0

si_x = 800
si_y = 800
si_z = 33

si_x1 = si_x + 1
si_y1 = si_y + 1
si_z1 = si_z + 1

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


dx1.astype(binprec).tofile('dx.box')
dy1.astype(binprec).tofile('dy.box')
dz1.astype(binprec).tofile('dz.box')


# ==== physical parameters
fMin = 3.78e-05
fMax = 1.3e-4

fmid = 0.5*(fMin + fMax)
  
beta = (fMax-fMin)/Ly
ff   = np.linspace(fMin,fMax,si_y)

print('f_south = {0}; beta = {1}'.format(fMin,beta) )


#%==================== LAND ===================================

landh  = np.zeros((si_y,si_x));

landh = -Lz + landh

# walls
landh[:,0] = 0.0
landh[-1,:] = 0.0
landh.astype(binprec).tofile('topog.box')

#%=============== Surface forcing ===================================
# -- temperature --
sst = np.zeros((si_y,si_x));

TS = 20.0
TN = 0.0

sst = (TN-TS)*yg/Ly + TS

#thetaClimFile
sst.astype(binprec).tofile('sstclim.box')

# -- salinity --
empmr = np.zeros((si_y,si_x));

F0 = 1e-7 # m/s ~ 10mm/day
empmr = F0*np.sin(2*np.pi*yg/Ly)
empmr.astype(binprec).tofile('empmr.box')


# -- wind -- 
windx = np.zeros((si_y,si_x));

tauW = 0.2
windx = -tauW*np.sin(2*np.pi*yg/Ly )

#windx = windx*ff.reshape(si_y,1)/fMin
windx.astype(binprec).tofile('windx.box')

#=============== Initial conditions ===================================

uinit = np.zeros((si_z,si_y,si_x));
vinit = np.zeros((si_z,si_y,si_x));
tinit = np.zeros((si_z,si_y,si_x))
sinit = 35 + np.zeros((si_z,si_y,si_x))
einit = np.zeros((si_y,si_x));

uinit.astype(binprec).tofile('uinit.box')
vinit.astype(binprec).tofile('vinit.box')
tinit.astype(binprec).tofile('tinit.box')
sinit.astype(binprec).tofile('sinit.box')
einit.astype(binprec).tofile('einit.box')

#% ================ floats =============================                      

nfl0 = 22
nfl = nfl0**2
xfl = np.linspace(1.5*dx,Lx-0.5*dx,nfl0) # no float in topography
yfl = np.linspace(0.5*dx,Ly-1.5*dx,nfl0)
xxfl,yyfl = np.meshgrid(xfl,yfl) 

var_fl = np.zeros((nfl+1,9))


# float id
var_fl[:,0] = np.linspace(0,nfl,nfl+1)
# tstart
var_fl[:,1] = 0.0
#xpart
var_fl[1:,2] = xxfl.flatten() # np.linspace(0,Lx,nfl)
#ypart
var_fl[1:,3] = yyfl.flatten() #np.linspace(0,Lx,nfl)
#kpart
var_fl[:,4] = 0.0
#kfloat
var_fl[:,5] = 5.0
#iup
var_fl[:,6] = 0.0
#itop
var_fl[:,7] = 0.0
#tend
var_fl[:,8] = -1


#first line
var_fl[0,0] = nfl*1.0
var_fl[0,1] = -1 
var_fl[0,5] = nfl*1.0
var_fl[0,8] = -1

var_fl.astype(binprec).tofile('flinit.box')
