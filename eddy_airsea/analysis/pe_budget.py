#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import MITgcmutils as mit


plt.ion()

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


dir0 = 'tmp_energy7/'

filets = 'diag_ocnSnap*'
filepe = 'tracer_wb*'
fileave = 'diag_ocnTave*'

flag_grid = 0

alphat = 2e-4
betas = 7.4e-4

#%==================== LOAD FIELDS ===================================
RC    = mit.rdmds(dir0+'RC*')
RA    = mit.rdmds(dir0+'RA*')
DRF   = mit.rdmds(dir0+'DRF*')
hFacC = mit.rdmds(dir0+'hFacC*')

si_z,si_y,si_x = hFacC.shape

hFacC2 = np.where(hFacC < 1, np.NaN, 1)
RA = RA[None,:,:]

i = 1
iterst = mit.mds.scanforfiles(dir0 + filets)
itersp = mit.mds.scanforfiles(dir0 + filepe)

# t0 = mit.rdmds(dir0 + filets,iterst[i],rec=0)
# t1 = mit.rdmds(dir0 + filets,iterst[i+1],rec=0)

# s0 = mit.rdmds(dir0 + filets,iterst[i],rec=1)
# s1 = mit.rdmds(dir0 + filets,iterst[i+1],rec=1)

#w0 = mit.rdmds(dir0 + filew,iterst[i],rec=0)
#w1 = mit.rdmds(dir0 + filew,iterst[i+1],rec=0)

wav = mit.rdmds(dir0 + fileave,itersp[i],rec=4) 

dtdt = mit.rdmds(dir0 + filepe,itersp[i],rec=0)
dsdt = mit.rdmds(dir0 + filepe,itersp[i],rec=1)

advrt = mit.rdmds(dir0 + filepe,itersp[i],rec=2)
advxt = mit.rdmds(dir0 + filepe,itersp[i],rec=3)
advyt = mit.rdmds(dir0 + filepe,itersp[i],rec=4)

advrs = mit.rdmds(dir0 + filepe,itersp[i],rec=5)
advxs = mit.rdmds(dir0 + filepe,itersp[i],rec=6)
advys = mit.rdmds(dir0 + filepe,itersp[i],rec=7)

wb2 = mit.rdmds(dir0 + filepe,itersp[i],rec=8)
wb = mit.rdmds(dir0 + filepe,itersp[i],rec=9)

dtdt = dtdt/86400
dsdt = dsdt/86400

# t0 = np.where(t0 == 0,np.NaN,t0)
# t1 = np.where(t1 == 0,np.NaN,t1)

ix = np.int(si_x/2)

advrt = np.append(advrt,advrt[None,0,:,:],axis=0)
advrs = np.append(advrs,advrs[None,0,:,:],axis=0)

advyt = np.append(advyt,advyt[:,None,0,:],axis=1)
advys = np.append(advys,advys[:,None,0,:],axis=1)

advxt = np.append(advxt,advxt[:,:,None,0],axis=2)
advxs = np.append(advxs,advxs[:,:,None,0],axis=2)

adv_at = -( advrt[:-1,:,:] - advrt[1:,:,:] \
  + advyt[:,1:,:] - advyt[:,:-1,:]      \
  + advxt[:,:,1:] - advxt[:,:,:-1])     \
  /(RA*DRF)

adv_as = -( advrs[:-1,:,:] - advrs[1:,:,:] \
  + advys[:,1:,:] - advys[:,:-1,:]      \
  + advxs[:,:,1:] - advxs[:,:,:-1])     \
  /(RA*DRF)

def comp_b (temp,salt):
  return alphat*temp + betas*salt

# b0 = comp_b(t0,s0)
# b1 = comp_b(t1,s1)

# pe0 = RC*b0
# pe1 = RC*b1

dbdt = comp_b (dtdt,dsdt)
advb = comp_b(adv_at, adv_as)
dpedt = RC*(dbdt-advb)


def yzplot(psi,*args, **kwargs):
  
  psi = np.where(np.isnan(psi),0.,psi)
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

# psi = (pe1-pe0)/4500
# psi = psi[:,:,ix]
# yzplot(psi,title=r"tottend (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-5)

# psi2 = dpedt[:,:,ix]
# yzplot(psi2,title=r"tottend (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-5)

# psi3 = (psi - psi2)/psi2
# yzplot(psi3,title=r"tottend (m\,s$^{-2}$)",fgrid=flag_grid,vmax = 1e-3)


psi4 = adv_a[:,:,ix]
yzplot(psi4,vmax = 1e-4)

psi5 = dtdt[:,:,ix]
yzplot(psi5,vmax = 1e-4)

psi6 = psi5 - psi4
yzplot(psi6,vmax = 1e-4)
