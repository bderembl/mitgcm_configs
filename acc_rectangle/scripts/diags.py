#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit

plt.ion()
#plt.switch_backend('agg')

#Sverdrup
Sv = 1e6

iexp = ['1']

dir0 = '../run_06_p'
dir1 = dir0 + iexp[0] + '/mnc_test_*/'
file_s = 'state*.nc'
file_m = 'mean*.nc'
file_g = 'grid*.nc'

dir_a = '../run10/'
file_a = 'average.nc'

ntmax = 1000

fg = mit.mnc_files(dir1 + file_g)
hfacc = fg.variables['HFacC'][:,:,:].copy()
hfacw = fg.variables['HFacW'][:,:,:].copy()
hfacs = fg.variables['HFacS'][:,:,:].copy()

xx = fg.variables['X'][:].copy()
yy = fg.variables['Y'][:].copy()
xp1 = fg.variables['Xp1'][:].copy()
yp1 = fg.variables['Yp1'][:].copy()
dxc = fg.variables['dxC'][:].copy()
dyc = fg.variables['dyC'][:].copy()
dxf = fg.variables['dxF'][:].copy()
dyf = fg.variables['dyF'][:].copy()
drf = fg.variables['drF'][:].copy()
zz = -fg.variables['Z'][:].copy()
zl = -fg.variables['Zl'][:].copy()

fg.close()

iz500 = np.argmin(np.abs(zz-500))
iz1000 = np.argmin(np.abs(zz-1000))

Tacc = np.zeros(ntmax)
Tave = np.zeros(ntmax)
Tave500 = np.zeros(ntmax)
Tave1000 = np.zeros(ntmax)

vol = (hfacc*dxf*dyf*drf[:,None,None]).sum()
vol500 = (hfacc*dxf*dyf*drf[:,None,None])[:iz500].sum()
vol1000 = (hfacc*dxf*dyf*drf[:,None,None])[:iz1000].sum()

nit = 0

nexp = len(iexp)
for  ie in range(0,nexp):

  dir1 = dir0 + iexp[ie] + '/mnc_test_*/'
  f1 = mit.mnc_files(dir1 + file_s)
  
  it = f1.variables['iter'][:].copy()
  si_t = len(it)
  
  
  it0 = 0
  # remove last point except for last exp
  if ie < nexp - 1:
    it1 = si_t-1
  else:
    it1 = si_t

  for it in range(it0,it1):
  
    uvel  = f1.variables['U' ][it,:,:,:].copy()
    vvel  = f1.variables['V' ][it,:,:,:].copy()
    temp  = f1.variables['Temp' ][it,:,:,:].copy()

    Tacc[nit] = (uvel[:,:,0]*hfacw[:,:,0]*dxf[:,0]*drf[:,None]).sum()/Sv
    Tave[nit] = (temp*hfacc*dxf*dyf*drf[:,None,None]).sum()/vol
    Tave500[nit] = (temp*hfacc*dxf*dyf*drf[:,None,None])[:iz500].sum()/vol500
    Tave1000[nit] = (temp*hfacc*dxf*dyf*drf[:,None,None])[:iz500].sum()/vol1000


    nit += 1

  f1.close()
  

Tacc = Tacc[:nit]
Tave = Tave[:nit]
Tave500 = Tave500[:nit]
Tave1000 = Tave1000[:nit]


# ugly vorticity snapshot
dx = xx[1] - xx[0]
dudy = np.diff(uvel[0],1,0)
dvdx = np.diff(vvel[0],1,1)
vort = (dvdx[:-2,:] - dudy[:,:-2])/dx

#plt.imshow(vort/1e-5, origin='lower', cmap=plt.cm.seismic, vmin = -5, vmax = 5)                    
#plt.colorbar()                                                                                     

