#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import MITgcmutils as mit


plt.ion()

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


dir0 = '../run/'
dir_grid = dir0 + 'mnc_test_0001/'

file1 = 'diagU*'
file2 = 'diagV*'
file3 = 'diagKEU*'
file4 = 'diagKEV*'
file5 = 'diagKEs*'
fileg = 'grid.t001.nc'

flag_uv = 1 # 1: u , 2: v, 3:KEu, 4:KEv, 5:KEE
flag_grid = 1

#%==================== LENGTH SCALES ===================================

L0 = 1000
u0z = 1e-3     # vertical shear (du/dz) s-1
U0 = 1.0
rho0 = 1e3
gg = 9.81      # gravity
rho0z = 9.2e-4   # stratification (d rho /dz)
alphaK = 2.0e-4
H0 = 99.5 # m not exactly 100m to avoid discretization issue 

nu = 1.0
J = gg*rho0z/(rho0*u0z**2)
S = H0/L0
delta = (nu/u0z/L0**2)**(1./3)
adimfac = delta*np.sqrt(J)*S**2/2



#%==================== LOAD FIELDS ===================================

# load grid
if flag_grid:
  gridm = mit.rdmnc(dir_grid + fileg)
  XC    = gridm['XC'   ]
  YC    = gridm['YC'   ]
  XG    = gridm['XG'   ]
  YG    = gridm['YG'   ]
  DXC   = gridm['dxC'  ]
  DYC   = gridm['dyC'  ]
  hFacC = gridm['HFacC']
  hFacS = gridm['HFacS']
  hFacW = gridm['HFacW']
  RAS   = gridm['rAs'  ]
  RAW   = gridm['rAw'  ]
  RAC   = gridm['rA'  ]
  RAZ   = gridm['rAz'  ]
  RC    = gridm['RC'   ]
  RF    = gridm['RF'   ]
  DRC   = gridm['drC'  ]
  DRF   = gridm['drF'  ]
  Depth = gridm['Depth']


if flag_uv == 1:
  filer = file1
elif flag_uv == 2:
  filer = file2
elif flag_uv == 3:
  filer = file3
elif flag_uv == 4:
  filer = file4
elif flag_uv == 5:
  filer = file5

i = -1
iters1 = mit.mds.scanforfiles(dir0 + filer)

utot   = mit.rdmds(dir0 + filer,iters1[i],rec=0)
uadv   = mit.rdmds(dir0 + filer,iters1[i],rec=1)
upress = mit.rdmds(dir0 + filer,iters1[i],rec=2)
u_eta  = mit.rdmds(dir0 + filer,iters1[i],rec=3)
udissh = mit.rdmds(dir0 + filer,iters1[i],rec=4)
udissv = mit.rdmds(dir0 + filer,iters1[i],rec=5)
uext   = mit.rdmds(dir0 + filer,iters1[i],rec=6)
u_ab   = mit.rdmds(dir0 + filer,iters1[i],rec=7)

utot = utot/86400

si_z,si_y,si_x = utot.shape
ix = np.int(si_x/2)


def yzplot(psi,*args, **kwargs):
  
  vmax = np.max(np.abs((psi)))
  vmax = kwargs.get('vmax', vmax)
  vmin = -vmax
  psi = np.where(psi<vmin,vmin,psi)
  psi = np.where(psi>vmax,vmax,psi)
  
  title = kwargs.get('title',None)

  fgrid = kwargs.get('fgrid', 0)

  if fgrid:
    xx = XC[0,:]*1e-3
    yy = RC[:]
  else:
    si_y,si_x = psi.shape
    xx = np.arange(si_x)
    yy = np.arange(si_y)
    
  plt.figure()
  plt.contourf(xx,yy,psi,100,cmap=plt.cm.seismic,vmin=vmin,vmax=vmax,extend='both')
  plt.colorbar(format='%.0e')
  plt.title(title)
  if fgrid:
    plt.xlabel('x (km)')
    plt.ylabel('z (m)')

# psi = utot[:,0,:]
# yzplot(psi,title=r"tottend (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-6)

psi = uadv + upress + udissv + udissh + u_eta + u_ab + uext
psi = psi[:,0,:]
yzplot(psi,title=r"sum (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-4)

psi = uext
psi = psi[:,0,:]
yzplot(psi,title=r"sum (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-4)



# error = np.abs(uadv) + np.abs(upress) + np.abs(udissv) + np.abs(udissh) + np.abs(u_eta) + np.abs(u_ab) + np.abs(uext)
# error2 = error[:,0,:]/np.abs(utot[:,0,:])

# psi3 = (psi - utot[:,0,:])/utot[:,0,:]/error2
# yzplot(psi3,title=r"relative error (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-6)


