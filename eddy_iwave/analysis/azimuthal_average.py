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

def interp_weights(xyz, uvw):
    naux,d = xyz.shape
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def interpolate(values, vtx, wts, fill_value=np.nan):
    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    ret[np.any(wts < 0, axis=1)] = fill_value
    return ret

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


dir0 = '/run/media/bderembl/girtab/eddy-iwave/run13/'
#dir0 = '/home/bderembl/work/MITgcm/mitgcm_configs/eddy_iwave/run/'

file1 = 'diagU*'
file2 = 'diagV*'
file3 = 'diagSurf*'
file4 = 'U*'
file5 = 'V*'
file6 = 'W*'

#%==================== LOAD GRID ===================================

nml = f90nml.read(dir0+'data')
nmldiag = f90nml.read(dir0+'data.diagnostics')

# load grid
XC    = mit.rdmds(dir0+'XC*')
YC    = mit.rdmds(dir0+'YC*')
XG    = mit.rdmds(dir0+'XG*')
YG    = mit.rdmds(dir0+'YG*')
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
f0 = nml['parm01']['f0']

hFacC2 = np.where(hFacC != 1, np.nan,1.0)
hFacS2 = np.where(hFacS != 1, np.nan,1.0)
hFacW2 = np.where(hFacW != 1, np.nan,1.0)

iz = np.argmin(np.abs(RC+np.min(Depth)))

si_z,si_y,si_x = hFacC.shape

Lx = XC[-1,-1] + XC[0,0]
Ly = YC[-1,-1] + YC[0,0]

dx = 2*XC[0,0]
dy = 2*YC[0,0]

xy_g = np.vstack((XG.flatten(), YG.flatten())).T
xy_u = np.vstack((XG.flatten(), YC.flatten())).T
xy_v = np.vstack((XC.flatten(), YG.flatten())).T
xy_c = np.vstack((XC.flatten(), YC.flatten())).T


iters1 = mit.mds.scanforfiles(dir0 + file1)
iters4 = mit.mds.scanforfiles(dir0 + file4)

# ==== eddy parameters (cf. mygendata) ========

x_c = Lx/2
y_c = Ly/2
R0 = 14e3
velmax = 0.1

rr    = np.linspace(0.0,0.5*Lx, np.int(0.5*si_x)+1)
theta = np.linspace(0.0,2*np.pi, np.int(np.pi*Lx/dx))
theta = theta[:-1]
rg,tg = np.meshgrid(rr,theta) 
si_t,si_r = rg.shape
dr = rr[1] - rr[0]
rr2 = rr[:-1] + rr[1:]


x_rt = rg*np.cos(tg) + x_c
y_rt = rg*np.sin(tg) + y_c
xy_rt = np.vstack((x_rt.flatten(), y_rt.flatten())).T

vtx_g, wts_g = interp_weights(xy_g, xy_rt)
vtx_u, wts_u = interp_weights(xy_u, xy_rt)
vtx_v, wts_v = interp_weights(xy_v, xy_rt)
vtx_c, wts_c = interp_weights(xy_c, xy_rt)

# grid at U,V,T points
rad_gg = np.sqrt((XG-x_c)**2 + (YG-y_c)**2)
rad_cc = np.sqrt((XC-x_c)**2 + (YC-y_c)**2)
rad_gu = np.sqrt((XG-x_c)**2 + (YC-y_c)**2)
rad_gv = np.sqrt((XC-x_c)**2 + (YG-y_c)**2)

theta_gg = np.arctan2(YG-y_c,XG-x_c)
theta_cc = np.arctan2(YC-y_c,XC-x_c)
theta_gu = np.arctan2(YC-y_c,XG-x_c)
theta_gv = np.arctan2(YG-y_c,XC-x_c)


# vortex
def vel_rankine(rr):
  v = -velmax*np.tanh(rr/R0)/(np.cosh(rr/R0))**2/(np.tanh(1.0)/(np.cosh(1.0))**2)
  v = np.where(rr == 0, 0.0,v)
  return v

#%==================== LOAD FIELDS ===================================


i = 1
udissv = mit.rdmds(dir0 + file1,iters1[i],rec=1)
vdissv = mit.rdmds(dir0 + file2,iters1[i],rec=1)

uvel = mit.rdmds(dir0 + file4,iters4[i])
vvel = mit.rdmds(dir0 + file5,iters4[i])
wvel = mit.rdmds(dir0 + file6,iters4[i])

#uvel0 = mit.rdmds(dir0 + file4,iters4[i-1])
#vvel0 = mit.rdmds(dir0 + file5,iters4[i-1])

ur_me = np.zeros((si_z,si_r))
ut_me = np.zeros((si_z,si_r))
ur_me1 = np.zeros((si_z,si_r))
ut_me1 = np.zeros((si_z,si_r))
w_me1 = np.zeros((si_z,si_r))

#ut0_me = np.zeros((si_z,si_r))


urdissv_me = np.zeros((si_z,si_r))
utdissv_me = np.zeros((si_z,si_r))

stress1 = np.zeros((si_z+1,si_r))
stress2 = np.zeros((si_z,si_r))
stress3 = np.zeros((si_z+1,si_r))

# set topography points to nans
uvel = uvel*hFacW2
vvel = vvel*hFacS2
wvel = wvel*hFacC2


for k in range(0,si_z-1):
  udissv[k,:,:] = (udissv[k+1,:,:] - udissv[k,:,:])/(RAW*DRF[k]*hFacW[k,:,:])
  vdissv[k,:,:] = (vdissv[k+1,:,:] - vdissv[k,:,:])/(RAS*DRF[k]*hFacS[k,:,:])
    
udissv[si_z-1,:,:] = 0.0
vdissv[si_z-1,:,:] = 0.0
    
for k in range(0,si_z-1):
  uvel_pol = interpolate(uvel[k,:,:], vtx_u, wts_u).reshape((si_t,si_r))
  vvel_pol = interpolate(vvel[k,:,:], vtx_v, wts_v).reshape((si_t,si_r))
  uvel_pol1 = interpolate(uvel[k+1,:,:], vtx_u, wts_u).reshape((si_t,si_r))
  vvel_pol1 = interpolate(vvel[k+1,:,:], vtx_v, wts_v).reshape((si_t,si_r))

  wvel_pol1 = interpolate(wvel[k+1,:,:], vtx_c, wts_c).reshape((si_t,si_r))

  udissv_pol = interpolate(udissv[k,:,:], vtx_u, wts_u).reshape((si_t,si_r))
  vdissv_pol = interpolate(vdissv[k,:,:], vtx_v, wts_v).reshape((si_t,si_r))

  # u and v at vertical cell face
  uvel_pol1 = 0.5*(uvel_pol + uvel_pol1)
  vvel_pol1 = 0.5*(vvel_pol + vvel_pol1)


  ur =  np.cos(tg)*uvel_pol + np.sin(tg)*vvel_pol 
  ut = -np.sin(tg)*uvel_pol + np.cos(tg)*vvel_pol 

  ur1 =  np.cos(tg)*uvel_pol1 + np.sin(tg)*vvel_pol1 
  ut1 = -np.sin(tg)*uvel_pol1 + np.cos(tg)*vvel_pol1 

  urdissv =  np.cos(tg)*udissv_pol + np.sin(tg)*vdissv_pol 
  utdissv = -np.sin(tg)*udissv_pol + np.cos(tg)*vdissv_pol 


  ur_me[k,:] = np.nanmean(ur,axis=0)
  ut_me[k,:] = np.nanmean(ut,axis=0)
  ur_me1[k,:] = np.nanmean(ur1,axis=0)
  ut_me1[k,:] = np.nanmean(ut1,axis=0)
  w_me1 [k,:] = np.nanmean(wvel_pol1,axis=0)
  urdissv_me[k,:] = np.nanmean(urdissv,axis=0)
  utdissv_me[k,:] = np.nanmean(utdissv,axis=0)


  # uvel_pol = interpolate(uvel0[k,:,:], vtx_u, wts_u).reshape((si_t,si_r))
  # vvel_pol = interpolate(vvel0[k,:,:], vtx_v, wts_v).reshape((si_t,si_r))
  # ut0 = -np.sin(tg)*uvel_pol + np.cos(tg)*vvel_pol 
  # ut0_me[k,:] = np.nanmean(ut0,axis=0)


  stress1[k+1,:] = -np.nanmean((ut1 - ut_me1[k,:])*(wvel_pol1 - w_me1[k,:]),axis=0)
  stress2[k,:] = -np.nanmean(rr.reshape((1,si_r))*(ut - ut_me[k,:])*(ur - ur_me[k,:]),axis=0)

# minus DRF because diff done downward
stressdiv1 = np.diff(stress1,1,0)/(-DRF[:,0,:])
stressdiv2 = 1/rr2.reshape((1,si_r-1))*np.diff(stress2,1,1)/dr

stressdiv = stressdiv1[:,1:] + stressdiv2

dutdz = np.diff(ut_me,1,0)/(-DRF[:-1,0,:])


#================ Plot part ================

def rzplot(psi,*args, **kwargs):
  
  vmax = np.max(np.abs((psi)))
  vmax = kwargs.get('vmax', vmax)
  vmin = -vmax
  psi = np.where(psi<vmin,vmin,psi)
  psi = np.where(psi>vmax,vmax,psi)
  
  title = kwargs.get('title',None)

  plt.figure()
  plt.contourf(rr*1e-3,RC[:,0,0],psi,100,cmap=plt.cm.seismic,vmin=vmin,vmax=vmax,extend='both') 
  plt.colorbar(format='%.0e')
  plt.contour(rr*1e-3,RC[:,0,0],ut_me,np.linspace(-0.2,0.2,17),colors='k',linewidths=0.5)
  plt.xlabel('r (km)')
  plt.ylabel('z (m)')
  plt.title(title)


vmaxall = 3e-7
psi  = -f0*ur_me
rzplot(psi,title=r"$-fU_r$ (m\,s$^{-2}$)",vmax=vmaxall)
plt.savefig('ucori.png',bbox_inches='tight')

# psi = rr*ut_me + 0.5*f0*rr**2
# rzplot(psi,title=r"$\lambda$ (m$^2$\,s$^{-1}$)")

# #psi = (ut_me-vel_rankine(rr))/(iters4[1]*dt)
# psi = (ut_me-ut0_me)/((iters4[1]-iters4[0])*dt)
# rzplot(psi,title=r"$du_\theta/dt$ (m\,s$^{-2}$)",vmax=vmaxall)

psi = stressdiv1
rzplot(psi,title=r"$\partial \overline{u'_\theta w'}/\partial z$ (m\,s$^{-2}$)",vmax=vmaxall)
plt.savefig('dupwpdz.png',bbox_inches='tight')

psi = utdissv_me
rzplot(psi,title=r"$\nu d^2 u_\theta/dz^2$ (m\,s$^{-2}$)",vmax=vmaxall)
plt.savefig('uvisc.png',bbox_inches='tight')

# vmin = -3e-7
# vmax = -vmin
# psi = stressdiv1[:iz-1,:]
# psi = np.where(psi<vmin,vmin,psi)
# psi = np.where(psi>vmax,vmax,psi)


# plt.figure()
# plt.contourf(rr*1e-3,RC[:iz-1,0,0],psi,100,cmap=plt.cm.seismic,vmin=vmin,vmax=vmax)
# plt.colorbar(format='%.0e',label='m/s2')
# plt.contour(rr*1e-3,RC[:,0,0],ut_me,5,colors='k',linewidths=0.5)
# plt.xlabel('r (km)')
# plt.ylabel('z (m)')
# plt.title(r"$\partial \overline{u'_\theta w'}/\partial z$")

# vmin = -3e-7
# vmax = -vmin
# #psi = (ut_me[:iz-1,:]-vel_rankine(rr))/(iters1[1]*dt)
# #psi = -(uvel[:iz-1,499:,500]-uvel0[:iz-1,499:,500])/(iters1[1]*dt)
# #psi = -(udiss[:iz-1,499:,500])



# plt.figure()
# plt.contourf(rr*1e-3,RC[:iz-1,0,0],psi,100,cmap=plt.cm.seismic,vmin=vmin,vmax=vmax,extend='both')
# plt.colorbar(format='%.0e',label='m/s2')
# plt.contour(rr*1e-3,RC[:,0,0],ut_me,5,colors='k',linewidths=0.5)
# plt.xlabel('r (km)')
# plt.ylabel('z (m)')
# plt.title(r"$\nu d^2 u_\theta/dz^2$")

# # fit

# # itaudrag = stressdiv1[:iz-1,:]/ut_me[:iz-1,:]
# # # mean over the eddy
# # itaudrag_me = itaudrag[:,50:250].mean(axis=1)
# # #plt.contourf(rr*1e-3,RC[:iz-1,0,0],itaudrag,100,cmap=plt.cm.seismic)
# # #plt.colorbar(format='%.0e',label='1/s')
# # plt.figure(); 
# # plt.plot(RC[:iz-1,0,0],1/(itaudrag_me*86400))
# # # fit
# # plt.plot(RC[:iz-1,0,0],1/(-0.3*(np.exp((-RC[:iz-1,0,0]-3800)/100)-1.6e-5*RC[:iz-1,0,0])))

# d2udz2 = np.diff(np.diff(uvel,1,0),1,0)/DRF[1:-1]**2

# nu = 1e-2
# plt.figure()
# plt.pcolormesh(nu*d2udz2[:iz-1,499:,500],cmap=plt.cm.seismic,vmin=vmin,vmax=vmax)
# plt.colorbar(format='%.0e',label='m/s2')
