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

file1 = 'diagU*'
file2 = 'diagV*'
file3 = 'diagKEU*'
file4 = 'diagKEV*'
file5 = 'diagKEs*'

flag_uv = 5 # 1: u , 2: v, 3:KEu, 4:KEv, 5:KEE
flag_grid = 1

#%==================== LOAD FIELDS ===================================

# load grid
if flag_grid:
  XC    = mit.rdmds(dir0+'XC*')
  YC    = mit.rdmds(dir0+'YC*')
  XG    = mit.rdmds(dir0+'XG*')
  YG    = mit.rdmds(dir0+'YG*')
  DXC   = mit.rdmds(dir0+'DXC*')
  DYC   = mit.rdmds(dir0+'DYC*')
  hFacC = mit.rdmds(dir0+'hFacC*')
  hFacS = mit.rdmds(dir0+'hFacS*')
  hFacW = mit.rdmds(dir0+'hFacW*')
  RAS   = mit.rdmds(dir0+'RAS*')
  RAW   = mit.rdmds(dir0+'RAW*')
  RAC   = mit.rdmds(dir0+'RAC*')
  RAZ   = mit.rdmds(dir0+'RAZ*')
  RC    = mit.rdmds(dir0+'RC*')
  RF    = mit.rdmds(dir0+'RF*')
  DRC   = mit.rdmds(dir0+'DRC*')
  DRF   = mit.rdmds(dir0+'DRF*')
  Depth = mit.rdmds(dir0+'Depth*')


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

i = 4
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
    xx = YC[:,ix]*1e-3
    yy = RC[:,0,0]
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

psi = utot[:,:,ix]
yzplot(psi,title=r"tottend (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-6)

psi = uadv + upress + udissv + udissh + u_eta + u_ab + uext
psi = psi[:,:,ix]
yzplot(psi,title=r"sum (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-6)

error = np.abs(uadv) + np.abs(upress) + np.abs(udissv) + np.abs(udissh) + np.abs(u_eta) + np.abs(u_ab) + np.abs(uext)
error2 = error[:,:,ix]/np.abs(utot[:,:,ix])

psi3 = (psi - utot[:,:,ix])/utot[:,:,ix]/error2
yzplot(psi3,title=r"relative error (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-7)


