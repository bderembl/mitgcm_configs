#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob,os,re
import sys

plt.ion()
if len(sys.argv) > 1:
  N = int(sys.argv[1])
else:
  print("usage: python mygendata.py N")
  print("with N the order of the grid (size of the grid is 2^N)")
  sys.exit(1)

# ==== physical constants
gg = 9.8      # gravity [m/s^2]
alphaT = 2e-4 # thermal expansion coef [1/K]
rho0 = 1000.0 # reference density [kg/m^3]
Cp   = 4000.0 # heat capacity [J/kg/K]

# ==== physical parameters
fmid = 8e-05 # Coriolis parameter in the middle [1/s]
beta = 2.e-11 # beta effect [1/m/s]
tauW = 0#0.08   # wind stress magnitude [N/m2]
TN = 0.0     # northern boundary temperature relaxation [K]
TS = 30.0    # southern boundary temperature relaxation [K]
qwatt = 35   # surface forcing [W/m2/K]
bottomDragLinear = 0.0002 # bottom drag [m/s]

# domain size 
Lx = 5000.0e3 # zonal length [m]
Ly = 5000.0e3 # meridional length [m]
Lz = 4000.0   # depth [m]

# vertical grid stretching parameter
stretch_coef = 4

# topography
wtopo   = 200e3  # width of the topography [m]
topomax = 2000    # maxheight of the topography [m]

#scales for non dimensional numbers
u_qg = 0.1    # characteristic velocity [m/s]
l_qg = 50e3   # characteristic length scale [m]

# ==== derived quantities: 
# number of grid points in x,y,z direction
si_x = 2**N
si_y = 2**N
si_z = 20 + int(max(si_x/32,1))

# number of processes
npx = int(max(si_x/64,1))
npy = int(max(si_x/64,1))

#viscosity
Uv_ref = 1.0     #  lateral viscous velocity [m/s]
Lv_ref = Lx/si_x #  lateral viscous length   [m] 

nu4 = 1./12*Uv_ref*Lv_ref**3
nuz = 1e-3*min(32/si_x,1) + 1e-5

# time stepping
dtcfl = 0.1*Lv_ref/Uv_ref     # cfl 
dtcflvisc = Lv_ref**4/nu4/32  # cfl on viscosity
deltat = min(dtcfl,dtcflvisc)
deltat = int(min(1800, deltat)/60)*60
dumpfreq = 360*86400    # [s]
endtime = 100*360*86400 # [s]

#binary precision
binprec = '>f4'

#% ================== GRID =====================================

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

# vertical grid
xfz = np.linspace(0,1,1000)
yfz = np.sinh(stretch_coef*xfz)/np.sinh(stretch_coef)

zc,zf,dz1 = stretch(xfz,yfz,Lz,si_z)

# Coriolis
ff = fmid + beta*(yg-Ly/2)

dx1.astype(binprec).tofile('dx.box')
dy1.astype(binprec).tofile('dy.box')
dz1.astype(binprec).tofile('dz.box')

#%==================== TOPO ===================================

# sides
def shelf(x,d):
  return (1-np.exp(-x**2/(2*d**2)))

topo = -Lz + topomax*(1-shelf(xg,wtopo)*shelf(Lx-xg,wtopo)*shelf(yg,wtopo)*shelf(Ly-yg,wtopo))

# walls
topo[:,0] = 0.0
topo[-1,:] = 0.0

topo.astype(binprec).tofile('topog.box')

#%=============== Surface forcing ===================================

# -- temperature --
sst = TS + (TN-TS)*yg/Ly
#sst = TS + (TN-TS)*np.sin(np.pi*yg/2/Ly)
#sst = TS + (TN-TS)*np.sin(np.pi*yg/2.5/Ly)/np.sin(np.pi/2.5)

#thetaClimFile
sst.astype(binprec).tofile('sstclim.box')

# -- wind -- 
windx = -tauW*ff/fmid*np.sin((2*np.pi*(yg-0.5*dx)/(Ly-dx)))

windx.astype(binprec).tofile('windx.box')

#%=============== initial conditions ===================================

uvel  = np.zeros((si_z,si_y,si_x));
vvel  = np.zeros((si_z,si_y,si_x));
theta = np.zeros((si_z,si_y,si_x));

eta   = np.zeros((si_y,si_x));

#theta = theta + temp_i - TN
#theta = theta*(1-yg/Ly) + TN


uvel.astype(binprec).tofile('uinit.box')
vvel.astype(binprec).tofile('vinit.box')
theta.astype(binprec).tofile('tinit.box')

eta.astype(binprec).tofile('einit.box')

#================= print interesting quantities

# relax time scale 
tauThetaClimRelax = rho0*Cp*dz1[0]/qwatt 
# rossby number
Ro = u_qg/(l_qg*fmid)
#thickness of ekman layer
de = 2*bottomDragLinear/fmid
#ekman number
Ek = de/Lz
taubottom = dz1[-1]/bottomDragLinear

print("Nb grid points: {0}x{0}x{1}".format(si_x,si_z))
print("qwatt = {0} W/m2/K".format(qwatt))
print("taurelax_sst = {0} days".format(tauThetaClimRelax/86400))
print("taubottom = {0} days".format(taubottom/86400))


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
      
with open("data_0", "r") as sources:
    lines = sources.readlines()
with open("data", "w") as sources:
    for line in lines:
      line2 = re.sub(r'MYDELTAT', str(deltat), line)
      line2 = re.sub(r'MYENDTIME', str(endtime), line2)
      line2 = re.sub(r'MYTAUTHETACLIMRELAX', str(tauThetaClimRelax), line2)
      line2 = re.sub(r'MYDUMPFREQ', str(dumpfreq), line2)
      line2 = re.sub(r'MYVISCA4', str(nu4), line2)
      line2 = re.sub(r'MYVISCAZ', str(nuz), line2)
      line2 = re.sub(r'MYDIFFK4T', str(nu4), line2)
      line2 = re.sub(r'MYDIFFKZT', str(nuz), line2)
      line2 = re.sub(r'MYBETA', str(beta), line2)
      line2 = re.sub(r'MYF0', str(fmid - beta*Ly/2), line2)
      line2 = re.sub(r'MYBOTTOMDRAGLINEAR', str(bottomDragLinear), line2)
      sources.write(line2)
