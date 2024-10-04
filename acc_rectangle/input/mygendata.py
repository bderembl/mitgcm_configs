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
rhoa = 1.2    # atmos reference density [kg/m^3]
Cp   = 4000.0 # heat capacity [J/kg/K]


# ==== physical parameters
fmid = -1.11e-4 # Coriolis parameter at southern edge [1/s]
beta = 1.47e-11 # beta effect [1/m/s]
tauW = 0.2   # wind stress magnitude [N/m2]
TN = 15.0     # northern boundary temperature relaxation [K]
TS = 0.0    # southern boundary temperature relaxation [K]
qwatt = 40   # surface forcing [W/m2/K]
bottomDragLinear = 1.1e-3 # bottom drag [m/s]
maxdt = 1800 # upper bound for time step [s]
Cdrag = 1e-3 #atmospheric drag coef (nondimensional)
useEXF = 1
useRelativeWind = 1

# domain size 
Lx = 9600.0e3 # zonal length [m]
Ly = 2400.0e3 # meridional length [m]
Lz = 3500.0   # depth [m]

# vertical grid stretching parameter
stretch_coef = 4

# topography
wtopo   = 250e3  # width of the topography [m]
topomax = 1500    # maxheight of the topography [m]

# rough topography
h_rough_rms = 250  # topography std  [m]
l1_rough = 100e3 # first wave length [m]
l2_rough = 300e3 # second wave length [m]


# ==== derived quantities: 
# number of grid points in x,y,z direction
si_x = 2**N
si_y = 2**(N-2)
si_z = 20 + int(max(si_x/32,1))

# number of processes
npx = int(max(si_x/64,1))
npy = int(max(si_y/64,1))

#viscosity
Uv_ref = 1.0     #  lateral viscous velocity [m/s]
Lv_ref = Lx/si_x #  lateral viscous length   [m] 

nu4 = 2e-2*Uv_ref*Lv_ref**3
nuz = 2e-4*min(32/si_x,1) + 1e-5

# time stepping
dtcfl = 0.1*Lv_ref/Uv_ref     # cfl 
dtcflvisc = Lv_ref**4/nu4/32  # cfl on viscosity
deltat = min(dtcfl,dtcflvisc)
deltat = int(min(maxdt, deltat)/60)*60
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

topo = -Lz + topomax*(1-shelf(xg-Lx/2,wtopo))

# create random topo
kx = np.fft.fftshift(np.fft.fftfreq(si_x,dx))
ky = np.fft.fftshift(np.fft.fftfreq(si_y,dy))
k,l = np.meshgrid(kx,ky)
K = np.sqrt(k**2 + l**2)

k_topo1 = 1./(2*l1_rough)
k_topo2 = 1./(2*l2_rough)
K_mask = np.where(K<k_topo2,0, K)
K_mask = np.where(K_mask>k_topo1,0, 1)

rough_topo_hat = np.random.rand(si_y,si_x)*np.exp(1j*2*np.pi*np.random.rand(si_y,si_x))
rough_topo_hat *= K_mask
rough_topo = np.fft.fft2(np.fft.ifftshift(rough_topo_hat)).real
rough_topo = h_rough_rms*rough_topo/np.std(2*rough_topo)

topo = topo + rough_topo + 2*h_rough_rms

# southern wall
topo[0,:] = 0.0

topo.astype(binprec).tofile('topog.box')

#%=============== Surface forcing ===================================

# -- temperature --
sst = TS + (TN-TS)*yg/Ly
#sst = TS + (TN-TS)*np.sin(np.pi*yg/2/Ly)
#sst = TS + (TN-TS)*np.sin(np.pi*yg/2.5/Ly)/np.sin(np.pi/2.5)

#thetaClimFile
sst.astype(binprec).tofile('sstclim.box')

# -- wind -- 
windx = tauW/2*(1. + np.cos((2*np.pi*(yg-0.5*dx - 0.5*Ly)/(Ly-dx))))
windvelx = np.sqrt(tauW/rhoa/Cdrag/2*(1. + np.cos((2*np.pi*(yg-0.5*dx - 0.5*Ly)/(Ly-dx)))))

windx.astype(binprec).tofile('windx.box')
windvelx.astype(binprec).tofile('windvelx.box')

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
#thickness of ekman layer
de = 2*bottomDragLinear/fmid
#ekman number
Ek = de/Lz
taubottom = dz1[-1]/bottomDragLinear

print("Nb grid points: {0}x{1}x{2}".format(si_x,si_y,si_z))
print("npx = {0}, npy = {1}".format(npx, npy))
print("qwatt = {0} W/m2/K".format(qwatt))
print("taurelax_sst = {0} days".format(tauThetaClimRelax/86400))
if (useRelativeWind):
  print("using Relative Wind: no bottom friction")
else:
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
      if useEXF:
        line2 = re.sub(r'MYTAUTHETACLIMRELAX', str(0.0), line2)
      else:
        line2 = re.sub(r'MYTAUTHETACLIMRELAX', str(tauThetaClimRelax), line2)
      line2 = re.sub(r'MYDUMPFREQ', str(dumpfreq), line2)
      line2 = re.sub(r'MYVISCA4', str(nu4), line2)
      line2 = re.sub(r'MYVISCAZ', str(nuz), line2)
      line2 = re.sub(r'MYDIFFK4T', str(nu4), line2)
      line2 = re.sub(r'MYDIFFKZT', str(nuz), line2)
      line2 = re.sub(r'MYBETA', str(beta), line2)
      line2 = re.sub(r'MYF0', str(fmid - beta*Ly/2), line2)
      if (useRelativeWind):
        line2 = re.sub(r'MYBOTTOMDRAGLINEAR', str(0.0), line2)
      else:
        line2 = re.sub(r'MYBOTTOMDRAGLINEAR', str(bottomDragLinear), line2)

      sources.write(line2)

with open("data.exf_0", "r") as sources:
    lines = sources.readlines()
with open("data.exf", "w") as sources:
    for line in lines:
      line2 = re.sub(r'MYCDRAG', str(Cdrag), line)
      if (useRelativeWind):
        line2 = re.sub(r'MYUSERELATIVEWIND', '.TRUE.', line2)
      else:
        line2 = re.sub(r'MYUSERELATIVEWIND', '.FALSE.', line2)
      line2 = re.sub(r'MYTAUTHETACLIMRELAX', str(tauThetaClimRelax), line2)
      sources.write(line2)

with open("data.pkg_0", "r") as sources:
    lines = sources.readlines()
with open("data.pkg", "w") as sources:
    for line in lines:
      if (useEXF):
        line2 = re.sub(r'MYUSEEXF', '.TRUE.', line)
      else:
        line2 = re.sub(r'MYUSEEXF', '.FALSE.', line)
      sources.write(line2)
