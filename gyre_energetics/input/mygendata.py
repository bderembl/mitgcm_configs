
#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

plt.ion()

binprec = '>f4'
flag_plot = 0


#% ================== GRID =====================================
rSphere = 6370.e3
deg2m = 2*np.pi*rSphere/360.0
gg = 9.8

si_x = 100
si_y = 100
si_z = 33


si_x1 = si_x + 1
si_y1 = si_y + 1
si_z1 = si_z + 1

Lx = 5000.0e3
Ly = 5000.0e3
Lz = 5000.0

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

# xf is % of grid points
xf = [0, 0.4, 0.6, 0.8, 0.9, 1]
#xf = [0, 0.5, 1]
# yf is % of thickness
yf = [0, 0.08, 0.14, 0.21, 0.4, 1]
#yf = [0, 0.5, 1]

hh = np.linspace(0,1,si_z1)
zf = Lz*np.interp(hh,xf,yf)

# smooth
nc = int(si_z/10)
if nc % 2 == 0:
  nc = nc + 1
zz2 = np.convolve(zf, np.ones((nc,))/nc, mode='valid')

zf[int((nc-1)/2):int(-(nc-1)/2)] = zz2

if flag_plot:
  plt.figure()
  plt.plot(hh,zz/Lz,'k')
  plt.plot(hh,hh,'k--')
  plt.plot(xf,yf,'.')
  plt.savefig(outputdir2 + 'vert_res.png')
  plt.close()

dz1 = np.diff(zf)

iz = np.argmin(np.abs(zf-500.0))

print ('dx= ', dx)
print ('min dz: ', np.min(dz1))
print ('max dz: ', np.max(dz1))
print ('nb layers above 500m:', iz, '/', si_z)

if np.sum(dz1 < 0) > 0:
  print ('you need you change the polynomial fit!')

zc = zf[0:-1] + 0.5*dz1

# # 1 layer configuration
# si_z = 1
# dz1 = np.zeros((si_z))
# dz1[0] = 4000.0

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
einit = np.zeros((si_y,si_x));

uinit.astype(binprec).tofile('uinit.box')
vinit.astype(binprec).tofile('vinit.box')
tinit.astype(binprec).tofile('tinit.box')
einit.astype(binprec).tofile('einit.box')
