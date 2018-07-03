#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import itertools

import MITgcmutils as mit
import f90nml


plt.ion()

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


dir0 = '/home/bderembl/work/MITgcm/mitgcm_configs/eddy_iwave/run/'

file1 = 'diagU*'
file2 = 'diagV*'
file3 = 'diagSurf*'


#%==================== LOAD FIELDS ===================================

nml = f90nml.read(dir0+'data')
nmldiag = f90nml.read(dir0+'data.diagnostics')

nmldiag = f90nml.read(dir0+'data.diagnostics')
#if nmldiag['diagnostics_list']['fields'][0][1] == 'VISrI_Um':

# load grid
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

dt = nml['parm03']['deltat']
dtdiag = nmldiag['diagnostics_list']['frequency'][0]

si_z,si_y,si_x = hFacC.shape

gravity = 9.81

hFacC2 = np.where(hFacC != 1, np.nan,1.0)
hFacS2 = np.where(hFacS != 1, np.nan,1.0)
hFacW2 = np.where(hFacW != 1, np.nan,1.0)

iters1 = mit.mds.scanforfiles(dir0 + file1)
iters2 = mit.mds.scanforfiles(dir0 + file2)

i = 0
udissh = mit.rdmds(dir0 + file1,iters1[i],rec=0)
udissv = mit.rdmds(dir0 + file1,iters1[i],rec=1)
uadv   = mit.rdmds(dir0 + file1,iters1[i],rec=2)
ucori  = mit.rdmds(dir0 + file1,iters1[i],rec=3)
uext   = mit.rdmds(dir0 + file1,iters1[i],rec=4)
upress = mit.rdmds(dir0 + file1,iters1[i],rec=5)
u_ab   = mit.rdmds(dir0 + file1,iters1[i],rec=6)
utot   = mit.rdmds(dir0 + file1,iters1[i],rec=7)

uvel0  = mit.rdmds(dir0 + file1,iters1[np.min([i-1,0])],rec=8)
uvel1  = mit.rdmds(dir0 + file1,iters1[i]  ,rec=8)

vvel = mit.rdmds(dir0 + file2,iters1[i]  ,rec=8)

psurf  = mit.rdmds(dir0 + file3,iters1[i],rec=0)

dpsdx = 0.*psurf
dpsdx[:,1:] = - gravity*(psurf[:,1:]-psurf[:,:-1])/DXC[:,:-1]

for k in range(0,si_z-1):
  udissv[k,:,:] = (udissv[k+1,:,:] - udissv[k,:,:])/(RAW*DRF[k]*hFacW[k,:,:])
#  vdissv[k,:,:] = (vdissv[k+1,:,:] - vdissv[k,:,:])/(RAS*DRF[k]*hFacS[k,:,:])
  
udissv[si_z-1,:,:] = 0.0
#vdissv[si_z-1,:,:] = 0.0


utot = utot/86400

ucori2 = 0.*ucori
ucori2[:,:-1,:] = 0.25*(vvel[:,1:,:] + vvel[:,:-1,:])
ucori2[:,:-1,1:] += 0.25*(vvel[:,1:,:-1] + vvel[:,:-1,:-1])

ucori2 = 1e-4*ucori2

ix = np.int(si_x/2)


def yzplot(psi,*args, **kwargs):
  
  vmax = np.max(np.abs((psi)))
  vmax = kwargs.get('vmax', vmax)
  vmin = -vmax
  psi = np.where(psi<vmin,vmin,psi)
  psi = np.where(psi>vmax,vmax,psi)
  
  title = kwargs.get('title',None)

  plt.figure()
  plt.contourf(YC[:,ix]*1e-3,RC[:,0,0],psi,100,cmap=plt.cm.seismic,vmin=vmin,vmax=vmax,extend='both')
  plt.colorbar(format='%.0e')
  plt.contour(YC[:,ix]*1e-3,RC[:,0,0],uvel1[:,:,ix],np.linspace(-0.2,0.2,17),colors='k',linewidths=0.5)
  plt.xlabel('r (km)')
  plt.ylabel('z (m)')
  plt.title(title)


#psi = uadv + upress + udissv + udissh + dpsdx + u_ab

psi = utot[:,:,ix]
yzplot(psi,title=r"utot (m\,s$^{-2}$)",vmax=1e-7)

psi = ucori2[:,:,ix]
yzplot(psi,title=r"ucori (m\,s$^{-2}$)",vmax=1e-7)


