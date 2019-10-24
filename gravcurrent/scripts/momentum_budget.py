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


dir0 = '../run/'

file1 = 'diagKEs*'
file2 = 'U*'
file3 = 'W*'
file4 = 'S*'


nml = f90nml.read(dir0+'data')
#nmldiag = f90nml.read(dir0+'data.diagnostics')

dt = nml['parm03']['dumpfreq']


#%==================== LOAD FIELDS ===================================

# load grid
XC    = mit.rdmds(dir0+'XC*').squeeze()
YC    = mit.rdmds(dir0+'YC*').squeeze()
XG    = mit.rdmds(dir0+'XG*').squeeze()
YG    = mit.rdmds(dir0+'YG*').squeeze()
DXC   = mit.rdmds(dir0+'DXC*').squeeze()
DYC   = mit.rdmds(dir0+'DYC*').squeeze()
hFacC = mit.rdmds(dir0+'hFacC*').squeeze()
hFacS = mit.rdmds(dir0+'hFacS*').squeeze()
hFacW = mit.rdmds(dir0+'hFacW*').squeeze()
RAS   = mit.rdmds(dir0+'RAS*').squeeze()
RAW   = mit.rdmds(dir0+'RAW*').squeeze()
RAC   = mit.rdmds(dir0+'RAC*').squeeze()
RAZ   = mit.rdmds(dir0+'RAZ*').squeeze()
RC    = mit.rdmds(dir0+'RC*').squeeze()
RF    = mit.rdmds(dir0+'RF*').squeeze()
DRC   = mit.rdmds(dir0+'DRC*').squeeze()
DRF   = mit.rdmds(dir0+'DRF*').squeeze()
Depth = mit.rdmds(dir0+'Depth*').squeeze()

si_z,si_x = hFacC.shape

gg = 9.81
sbeta = 7.4e-4  # psu^-1

hFacC2 = np.where(hFacC != 1, np.nan,1.0)
hFacS2 = np.where(hFacS != 1, np.nan,1.0)
hFacW2 = np.where(hFacW != 1, np.nan,1.0)

iters1 = mit.mds.scanforfiles(dir0 + file1)
iters2 = mit.mds.scanforfiles(dir0 + file2)

si_t = len(iters1)

s = mit.rdmds(dir0 + file4,iters2[0],rec=0).squeeze()
# parameters of the experiment
s1 = s[1,1]   # salty water
s0 = s[-1,-1] # fresh water
ipos = si_x - np.argmax(s[-1,::-1])
L0 = XC[ipos]
H0 = -RC[-1] - RC[0]
Lt = XC[-1] - XC[0]
dx = XC[1] - XC[0] # assume uniform grid
dz = RC[0] - RC[1] # assume uniform grid
gp = gg*sbeta*(s1-s0)
uref = np.sqrt(gp*H0)
tref = L0/uref

tt = dt*np.arange(0,si_t)

n0 = 0
n1 = si_t
uadv_s   = np.zeros(n1-n0)
upress_s = np.zeros(n1-n0)
upres2_s = np.zeros(n1-n0)
udissh_s = np.zeros(n1-n0)
udissv_s = np.zeros(n1-n0)
uext_s   = np.zeros(n1-n0)
u_ab_s   = np.zeros(n1-n0)
ke_s     = np.zeros(n1-n0)
dke_s    = np.zeros(n1-n0)
dkedt_s  = np.zeros(n1-n0)

def interp_z(w):
  w = np.pad(w,((0,1),(0,0)))
  w = 0.5*(w[:-1,:] + w[1:,:]) 
  return w

def interp_x(u):
  u = np.pad(u,((0,0),(0,1)))
  u = 0.5*(u[:,:-1] + u[:,1:]) 
  return u

for i in range(n0,n1):
  utot   = mit.rdmds(dir0 + file1,iters1[i],rec=0).squeeze()
  uadv   = mit.rdmds(dir0 + file1,iters1[i],rec=1).squeeze()
  upress = mit.rdmds(dir0 + file1,iters1[i],rec=2).squeeze()
  upres2 = mit.rdmds(dir0 + file1,iters1[i],rec=3).squeeze()
  udissh = mit.rdmds(dir0 + file1,iters1[i],rec=4).squeeze()
  udissv = mit.rdmds(dir0 + file1,iters1[i],rec=5).squeeze()
  uext   = mit.rdmds(dir0 + file1,iters1[i],rec=6).squeeze()
  u_ab   = mit.rdmds(dir0 + file1,iters1[i],rec=7).squeeze()
  
  u0 = mit.rdmds(dir0 + file2, iters1[i-1],rec=0).squeeze()
  u1 = mit.rdmds(dir0 + file2, iters1[i]  ,rec=0).squeeze()
  
  w0 = mit.rdmds(dir0 + file3, iters1[i-1],rec=0).squeeze()
  w1 = mit.rdmds(dir0 + file3, iters1[i]  ,rec=0).squeeze()
    
  dkedt = (interp_z(0.5*(w1**2 - w0**2)) + interp_x(0.5*(u1**2 - u0**2)))/dt
  ke = (interp_z(0.5*w1**2) + interp_x(0.5*u1**2))
  dke = uadv + upress + + upres2 + udissv + udissh + u_ab
  
  utot = utot/86400

  uadv_s  [i] = np.sum(uadv  )*dx*dz
  upress_s[i] = np.sum(upress)*dx*dz
  upres2_s[i] = np.sum(upres2)*dx*dz
  udissh_s[i] = np.sum(udissh)*dx*dz
  udissv_s[i] = np.sum(udissv)*dx*dz
  uext_s  [i] = np.sum(uext  )*dx*dz
  u_ab_s  [i] = np.sum(u_ab  )*dx*dz
  ke_s    [i] = np.sum(ke    )*dx*dz
  dkedt_s [i] = np.sum(dkedt )*dx*dz
  dke_s   [i] = np.sum(dke   )*dx*dz

# plt.figure()
# plt.plot(uadv_s  ,label='adv')
# plt.plot(upress_s,label='hd')
# plt.plot(upres2_s,label='nh')
# plt.plot(udissh_s,label='dish')
# plt.plot(udissv_s,label='disv')
# plt.plot(u_ab_s  ,label='ab')
# plt.legend()

plt.figure()
plt.plot(tt[n0:n1]/tref,ke_s/(uref**2*H0*L0),label='ke')
plt.plot(tt[n0:n1]/tref,np.cumsum(upress_s)*dt/(uref**2*H0*L0),label='hd')
plt.plot(tt[n0:n1]/tref,np.cumsum(upres2_s)*dt/(uref**2*H0*L0),label='nh')
plt.plot(tt[n0:n1]/tref,np.cumsum(udissh_s)*dt/(uref**2*H0*L0),label='dish')
plt.plot(tt[n0:n1]/tref,np.cumsum(udissv_s)*dt/(uref**2*H0*L0),label='disv')
#plt.plot(tt[n0:n1]/tref,np.cumsum(u_ab_s + uadv_s )*dt,label='num')
#plt.plot(np.cumsum(dke_s   )*dt,label='all')
#plt.plot(np.cumsum(dkedt_s )*dt,label='all2')
plt.legend()
plt.grid()

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


# plt.figure()
# plt.imshow(dkedt)    
# plt.colorbar()

# plt.figure()
# plt.imshow(dke)    
# plt.colorbar()

# plt.figure()
# plt.imshow(dke - dkedt)    
# plt.colorbar()


# plt.figure()
# plt.loglog(tt[n0:n1]/tref,np.cumsum(upress_s)*dt/(uref**2*H0*L0),label='hd')
# plt.loglog(tt[n0:n1]/tref,ke_s/(uref**2*H0*L0),label='ke')
# plt.grid()
