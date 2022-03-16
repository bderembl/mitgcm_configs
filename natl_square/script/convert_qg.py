#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit
import scipy.io.netcdf as netcdf
from scipy import interpolate
from scipy.ndimage.filters import gaussian_filter
import sys
import qgutils as qg
from scipy import ndimage as nd

import time

plt.ion()

if len(sys.argv) > 1:
  nl = int(sys.argv[1])
else:
  print("usage: python mygendata.py nl")
  print("with nl the number of qg layers")
  sys.exit(1)

# ==== physical constants
gg = 9.8      # gravity [m/s^2]
alphaT = 2e-4 # thermal expansion coef [1/K]
rho0 = 1000.0 # reference density [kg/m^3]
Cp   = 4000.0 # heat capacity [J/kg/K]

#scales for non dimensional numbers
u_qg = 0.1    # characteristic velocity [m/s]
l_qg = 50e3   # characteristic length scale [m]

# topography
wtopo   = 200e3  # width of the topography [m]
topomax = 2000    # maxheight of the topography [m]

smooth_rad = 50e3 # smoothing radius [m]
N2_min = 1e-7    # minimum N2 [1/s^2]

# vertical stretching coef for the QG layers
stretch_coef = 4

binprec = '>f4'

def fill(data, invalid=None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid') 
    by the value of the nearest valid data cell

    Input:
        data:    numpy array of any dimension
        invalid: a binary array of same shape as 'data'. True cells set where data
                 value should be replaced.
                 If None (default), use: invalid  = np.isnan(data)

    Output: 
        Return a filled array. 
    """
    #import numpy as np
    #import scipy.ndimage as nd

    if invalid is None: invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]

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




# average
dir0 = '../run10/'
dir1 = '../run10/'
file1 = 'average.nc'
file_g = 'grid.nc'

#f1 = netcdf.netcdf_file(dir1 + file1,'r')

if file1 == 'average.nc':
  f1 = netcdf.netcdf_file(dir1 + file1)
  fg = netcdf.netcdf_file(dir1 + file_g)
else:
  f1 = mit.mnc_files(dir1 + file1)
  fg = mit.mnc_files(dir1 + file_g)

hfacc = fg.variables['HFacC'][:,:,:].copy()
hfacw = fg.variables['HFacW'][:,:,:-1].copy()
hfacs = fg.variables['HFacS'][:,:-1,:].copy()
fcori = fg.variables['fCori'][:,:].copy()
xx = fg.variables['X'][:].copy()
yy = fg.variables['Y'][:].copy()
zz = -fg.variables['Z'][:].copy()
zl = -fg.variables['Zl'][:].copy()
fg.close()

it = f1.variables['iter'][:].copy()
si_t = len(it)

xg,yg = np.meshgrid(xx,yy) 


Lx = xx[-1] + xx[0]
Ly = yy[-1] + yy[0]
Lz = 2*zz[-1]-zl[-1]

si_x = len(xx)
si_y = len(yy)
si_z = len(zz)

dx = xx[1] - xx[0]
dy = yy[1] - yy[0]

zl = np.append(zl,Lz)
dzl = np.diff(zl)

uvel_me  = np.zeros((si_z,si_y,si_x))
vvel_me  = np.zeros((si_z,si_y,si_x))
theta_me = np.zeros((si_z,si_y,si_x))
#eta_me   = np.zeros((si_y,si_x))

nme = 0
it0 = si_t - 1
it1 = si_t
for it in range(it0,it1):

  uvel  = f1.variables['U'   ][it,:,:,:].copy()
  vvel  = f1.variables['V'   ][it,:,:,:].copy()
  theta = f1.variables['Temp'][it,:,:,:].copy()
#  eta   = f1.variables['Eta' ][it,:,:]

  uvel_me  += uvel[:si_z,:si_y,:si_x]
  vvel_me  += vvel[:si_z,:si_y,:si_x]
  theta_me += theta[:si_z,:si_y,:si_x]
#  eta_me   += eta[:si_y,:si_x]
  nme += 1


uvel_me  /= it1 - it0
vvel_me  /= it1 - it0
theta_me /= it1 - it0
#eta_me   /= it1 - it0

varu = uvel_me
varv = vvel_me
vart = theta_me
rhoa = -rho0*alphaT*theta_me 


# vertical grid
xfz = np.linspace(0,1,1000)
yfz = np.sinh(stretch_coef*xfz)/np.sinh(stretch_coef)
zc,zf,dh2 = stretch(xfz,yfz,Lz,nl)

si_z3 = len(dh2)
iz3 = np.zeros(si_z3+1,dtype='int')
dh3 = np.zeros(si_z3)
dh2c = dh2.cumsum()
dh2i = 0.5*(dh2[1:] + dh2[:-1])



z_interf = np.zeros(si_z3+1)
z_interf[1:] = dh2c

for nz in range(0,si_z3):
  iz3[nz+1] = np.argmin(np.abs(zl-dh2c[nz]))
  dh3[nz] = np.sum(dzl[iz3[nz]:iz3[nz+1]])

u_la = np.zeros((si_z3,si_y,si_x))
v_la = np.zeros((si_z3,si_y,si_x))
r_la = np.zeros((si_z3,si_y,si_x))

bz_la = np.zeros((si_z3,si_y,si_x))
with np.errstate(divide='ignore',invalid='ignore'):
  for nz in range(0,si_z3):
    u_la[nz,:,:] =  np.sum(hfacw[iz3[nz]:iz3[nz+1],:,:]*varu[iz3[nz]:iz3[nz+1],:,:]*dzl[iz3[nz]:iz3[nz+1],None,None],0)/(np.sum(hfacw[iz3[nz]:iz3[nz+1],:,:]*dzl[iz3[nz]:iz3[nz+1],None,None],0))
    v_la[nz,:,:] =  np.sum(hfacs[iz3[nz]:iz3[nz+1],:,:]*varv[iz3[nz]:iz3[nz+1],:,:]*dzl[iz3[nz]:iz3[nz+1],None,None],0)/(np.sum(hfacs[iz3[nz]:iz3[nz+1],:,:]*dzl[iz3[nz]:iz3[nz+1],None,None],0))
    r_la[nz,:,:] =  np.sum(hfacc[iz3[nz]:iz3[nz+1],:,:]*rhoa[iz3[nz]:iz3[nz+1],:,:]*dzl[iz3[nz]:iz3[nz+1],None,None],0)/(np.sum(hfacc[iz3[nz]:iz3[nz+1],:,:]*dzl[iz3[nz]:iz3[nz+1],None,None],0))

# adjust field in topography
u_la = np.where(np.isnan(u_la),0.,u_la)
v_la = np.where(np.isnan(v_la),0.,v_la)
   
for nz in range(0,si_z3):
  r_la[nz,:,:] = fill(r_la[nz,:,:])

N2 = gg/rho0*np.diff(r_la,1,0)/dh2i[:,None,None]
N2 = np.where(N2 < N2_min,N2_min,N2)

if smooth_rad > 0:
  print('Smoothing field')
  smooth_fac = smooth_rad/dx
  for nz in range(0,si_z3):
    u_la[nz,:,:] = gaussian_filter(u_la[nz,:,:],smooth_fac)
    v_la[nz,:,:] = gaussian_filter(v_la[nz,:,:],smooth_fac)
  for nz in range(0,si_z3-1):
    N2[nz,:,:] = gaussian_filter(N2[nz,:,:],smooth_fac)

print('Compute Rd')
#Rdqg = qg.comp_modes(dh3,N2,fcori,diag=True)


print('Solve stream function')

# using free slip BC
psi_ls = np.zeros((si_z3,si_y,si_x))
zeta = np.zeros((si_z3,si_x+1,si_x+1))
zeta[:,1:-1,1:-1] = (v_la[:,1:,1:] - v_la[:,1:,:-1])/dx - (u_la[:,1:,1:] - u_la[:,:-1,1:])/dy
zeta = 0.25*(zeta[:,:-1,:-1] + zeta[:,1:,1:] + zeta[:,:-1,1:] + zeta[:,1:,:-1])

for nz in range(0,si_z3):
  psi_ls[nz,:,:] = qg.solve_mg(zeta[nz,:,:],dx) 

U,V = qg.comp_vel(psi_ls, dx)

# non dimensional variables
Fr = u_qg/(np.sqrt(N2)*Lz)
psi_ls_a = psi_ls/(l_qg*u_qg)
#rda = Rdqg[1,:,:]/l_qg
dh_a = dh3/Lz

# topo
# sides
def shelf(x,d):
  return (1-np.exp(-x**2/(2*d**2)))

topo = topomax/Lz*(1-shelf(xg,wtopo)*shelf(Lx-xg,wtopo)*shelf(yg,wtopo)*shelf(Ly-yg,wtopo))


# Export to basilisk format
N = si_x
psi_ls_o = np.zeros((si_z3,N+1,N+1))
psi_ls_o[:,1:,1:] = psi_ls_a
psi_ls_o[:,0,:] = 0
psi_ls_o[:,:,0] = 0
psi_ls_o[:,0,0] = N
psi_ls_o = np.transpose(psi_ls_o,(0,2,1))

fileppg = dir0 + 'psipg_' + str(si_z3) +'l_N' + str(N) + '.bas'
psi_ls_o.astype('f4').tofile(fileppg)

print(f"Mean Fr = {np.mean(Fr,axis=(1,2))}")

Fr_o = np.zeros((si_z3,N+1,N+1))
Fr_o[:-1,1:,1:] = Fr
Fr_o[:,0,:] = 0
Fr_o[:,:,0] = 0
Fr_o[:,0,0] = N
Fr_o = np.transpose(Fr_o,(0,2,1))
fileFr = dir0 + 'frpg_' + str(si_z3) +'l_N' + str(N) + '.bas'
Fr_o.astype('f4').tofile(fileFr)

# # first deformation radius
# Rd_o = np.zeros((N+1,N+1))
# Rd_o[1:,1:] = rda
# Rd_o[0,:] = 0
# Rd_o[:,0] = 0
# Rd_o[0,0] = N
# Rd_o = np.transpose(Rd_o,(1,0))
# fileRd = dir0 + 'rdpg_' + str(si_z3) +'l_N' + str(N) + '.bas'
# Rd_o.astype('f4').tofile(fileRd)

# topography
topo_o = np.zeros((N+1,N+1))
topo_o[1:,1:] = topo
topo_o[0,:] = 0
topo_o[:,0] = 0
topo_o[0,0] = N
topo_o = np.transpose(topo_o,(1,0))
filetopo = dir0 + 'topo.bas'
topo_o.astype('f4').tofile(filetopo)


fileh = dir0 + 'dh_' + str(si_z3) +'l.bin'
dh_a.astype('f4').tofile(fileh)
print(f"dh = {dh_a}")

nx = int(si_x/3)
ny = int(si_y*2/3)

vmin = np.min(vart)
vmax = np.max(vart)

vcont = np.linspace(3,33,11) 

xkm = xx*1e-3
ykm = yy*1e-3

plt.figure()
CS = plt.contour(xkm,ykm,vart[0,:,:],vcont,colors='k')
plt.plot([xkm[nx], xkm[nx]],[ykm[0], ykm[-1]],'r',linewidth=2)
plt.plot([xkm[0], xkm[-1]],[ykm[ny], ykm[ny]],'r--',linewidth=2)
plt.clabel(CS)
plt.xlabel('x (km)')
plt.ylabel('y (km)')

idz = np.argmin(np.abs(zz-2000))

plt.figure()
ax = plt.subplot(2,1,1)
CS = plt.contourf(xkm,-zz[:],vart[:,ny,:],vcont, cmap=plt.cm.RdYlBu_r)
plt.contour(xkm,-zz[:idz],vart[:idz,ny,:],vcont,colors='w',linewidths=0.5)
plt.text(ykm[-18],-zz[idz-1],'EW section')
plt.xlabel("x (km)")
plt.ylabel("z (m)")
plt.xlim([0,5000])
#plt.clabel(CS)
ax = plt.subplot(2,1,2)
CS = plt.contourf(ykm,-zz[:],vart[:,:,nx],vcont, cmap=plt.cm.RdYlBu_r)
plt.contour(ykm,-zz[:idz],vart[:idz,:,nx],vcont,colors='w',linewidths=0.5)
plt.text(ykm[-18],-zz[idz-1],'NS section')
#plt.clabel(CS)
plt.xlabel("y (km)")
plt.ylabel("z (m)")
plt.xlim([0,5000])
plt.tight_layout()
plt.savefig('t_sec_pg.pdf')

# vcont = [5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0]
# plt.figure()
# plt.contourf(xkm,ykm,vart[0,:,:],25, cmap=plt.cm.RdYlBu_r)
# plt.colorbar()
# CSC = plt.contour(xkm,ykm,Rdqg[1,:,:]*1e-3,vcont,colors='k')
# plt.clabel(CSC)
# plt.xlabel('x (km)')
# plt.ylabel('y (km)')
# plt.savefig('sst_rd_pg.pdf')
# #plt.savefig('def_rad.png',bbox_inches='tight')

plt.figure()
plt.clf()
CS = plt.contour(xkm,ykm,psi_ls[0,:,:],10,colors ='k',linewidths=0.5)
ci = CS.levels[1] - CS.levels[0]
plt.text(3500,100,'ci:{0:.1e} m2/s'.format(ci*l_qg*u_qg))
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.xlim([0,5000])
plt.ylim([0,5000])
plt.savefig('ppg.pdf')
