#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
import MITgcmutils as mit
import qgutils as qg

plt.ion()

iexp = ['15','16']

n_coarse = 0

dir0 = '../run10_p'
dir1 = dir0 + iexp[0] + '/mnc_test_*/'
file_s = 'state*.nc'
file_m = 'mean*.nc'
file_g = 'grid*.nc'

dir_a = '../run10/'
file_a = 'average.nc'

fg = mit.mnc_files(dir1 + file_g)
hfacc = fg.variables['HFacC'][:,:,:].copy()
hfacw = fg.variables['HFacW'][:,:,:].copy()
hfacs = fg.variables['HFacS'][:,:,:].copy()

xx = fg.variables['X'][:].copy()
yy = fg.variables['Y'][:].copy()
xp1 = fg.variables['Xp1'][:].copy()
yp1 = fg.variables['Yp1'][:].copy()
zz = -fg.variables['Z'][:].copy()
zl = -fg.variables['Zl'][:].copy()

fg.close()


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




xg,yg = np.meshgrid(xx,yy) 

Lx = xp1[-1]
Ly = yp1[-1]
Lz = 2*zz[-1]-zl[-1]

si_x = len(xx)
si_y = len(yy)
si_z = len(zz)

dx = xx[1] - xx[0]
dy = yy[1] - yy[0]

dzl = np.diff(zl)

Nlev = int(np.log2(si_x))

si_xc = int(si_x/2**n_coarse)
si_yc = int(si_y/2**n_coarse)

f1 = netcdf.netcdf_file(dir_a + file_a)
uvel_me  = f1.variables['U'   ][0,:,:,:].copy()
vvel_me  = f1.variables['V'   ][0,:,:,:].copy()
f1.close()

if n_coarse > 0:
  uvel_me = qg.coarsen(uvel_me,n_coarse)
  vvel_me = qg.coarsen(vvel_me,n_coarse)


# QG vertical layers
#dh2 = np.array([300, 600, 1100, 2000])
#dh2 = dzl
# vertical grid
nl = si_z
# vertical grid stretching parameter
stretch_coef = 4
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


eke_me  = np.zeros((si_z3,si_yc,si_xc))

nme = 0

for  ie in range(0,len(iexp)):

  dir1 = dir0 + iexp[ie] + '/mnc_test_*/'
  f1 = mit.mnc_files(dir1 + file_s)
  
  it = f1.variables['iter'][:].copy()
  si_t = len(it)
  
  
  it0 = 0
  it1 = si_t-1
  for it in range(it0,it1):
  
    uvel  = f1.variables['U' ][it,:si_z,:si_y,:si_x].copy()
    vvel  = f1.variables['V' ][it,:si_z,:si_y,:si_x].copy()

    if n_coarse > 0:
      uvel = qg.coarsen(uvel,n_coarse)
      vvel = qg.coarsen(vvel,n_coarse)


    varu = uvel - uvel_me
    varv = vvel - vvel_me

    u_la = np.zeros((si_z3,si_yc,si_xc))
    v_la = np.zeros((si_z3,si_yc,si_xc))

    for nz in range(0,si_z3):
      u_la[nz,:,:] =  np.sum(varu[iz3[nz]:iz3[nz+1],:,:]*dzl[iz3[nz]:iz3[nz+1],None,None],0)/(np.sum(dzl[iz3[nz]:iz3[nz+1],None,None],0))
      v_la[nz,:,:] =  np.sum(varv[iz3[nz]:iz3[nz+1],:,:]*dzl[iz3[nz]:iz3[nz+1],None,None],0)/(np.sum(dzl[iz3[nz]:iz3[nz+1],None,None],0))

    eke = 0.5*(u_la**2 + v_la**2)
    eke = qg.wavelet(eke)

    eke_me += 0.5*(u_la**2 + v_la**2)
    nme += 1
  
  f1.close()
  
eke_me  /= nme

ke_me = 0.5*(uvel_me**2 + vvel_me**2)

plt.figure()
CS = plt.imshow(ke_me[0,:,:],origin='lower',vmax=0.1,cmap=plt.cm.hot,extent=(0,Lx*1e-3,0,Ly*1e-3))
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.colorbar()
plt.savefig('ke' + str(si_xc) + '.pdf')

plt.figure()
CS = plt.imshow(eke_me[0,:,:],origin='lower',vmax=0.1,cmap=plt.cm.hot,extent=(0,Lx*1e-3,0,Ly*1e-3))
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.colorbar()
plt.savefig('eke' + str(si_xc) + '.pdf')


eke_l = np.zeros((Nlev,si_y,si_x))

for il in range(0,Nlev):
  eke_l[il,:,:] = qg.wavelet_lowpass(eke[0,:,:],il)


plt.figure(); iw = 4
CS = plt.imshow(qg.refine(ekew[iw][0,:,:],Nlev-iw),origin='lower',vmax=0.1,cmap=plt.cm.hot,extent=(0,Lx*1e-3,0,Ly*1e-3))
