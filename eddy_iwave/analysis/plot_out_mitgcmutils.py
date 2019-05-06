#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit
import matplotlib

plt.ion()

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


dir0 = '/run/media/bderembl/girtab/eddy-iwave/run08/'

file1 = 'U*'
file2 = 'oceDiag*'


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


si_x = len(XC)
si_y = len(YC)
si_z = len(RC)

iters1 = mit.mds.scanforfiles(dir0 + file1)
iters2 = mit.mds.scanforfiles(dir0 + file2)

si_t1 = len(iters1)
si_t2 = len(iters2)

xmin = 1e10
xmax = -1e10
diss_me = np.zeros(si_z)
#diss_me2 = np.zeros(si_z)
nme = 0

#for i in iters2:
for i in [iters2[1]]:
  udiss = mit.rdmds(dir0 + file2,i,rec=1)
  uvel = mit.rdmds(dir0 + file2,i,rec=3)
  
  jumax,iumax = np.unravel_index(uvel[0,:,:].argmax(), uvel[0,:,:].shape)
  
  # if vertical dissipation
  for k in range(0,si_z-1):
    udiss[k,:,:] = (udiss[k+1,:,:] - udiss[k,:,:])/(RAW*DRF[k]*hFacW[k,:,:])

  diss = (udiss*uvel*hFacW*RAW).sum(axis=(1,2))

#  diss2 = udiss[:,400,500]*uvel[:,400,500]*RAW[400,500]

  diss_me = diss_me + diss  
#  diss_me2 = diss_me2 + diss2
  nme = nme + 1
  
diss_me = diss_me/nme/(hFacW*RAW).sum(axis=(1,2))
#diss_me2 = diss_me2/nme/(hFacW[:,400,500]*RAW[400,500])
xmin = np.nanmin([np.nanmin(diss_me),xmin])
xmax = np.nanmax([np.nanmax(diss_me),xmax])

plt.figure()
plt.plot([xmin, xmax],[-Depth.min(),-Depth.min()],'r--',linewidth=1)
plt.plot(diss_me,RC[:,0,0],'k',linewidth=1)
#plt.plot(diss_me2,RC[:,0,0],'k--',linewidth=1)
plt.xlabel(r'Dissipation (m$^4$s$^{-3}$)')
plt.ylabel(r'Height (m))')

#plt.savefig('dissip_wtopo2.pdf',bbox_inches='tight')


# uvel0 = mit.rdmds(dir0 + file1,0)
# plt.figure()
# plt.contourf(1e-3*YC[:,np.int(si_x/2)],RC[:,0,0],uvel0[:,:,np.int(si_x/2)],30,cmap=plt.cm.seismic)
# plt.colorbar(label=r'u (m\,s$^{-1}$)')
# plt.xlabel('y(km)')
# plt.ylabel('z(m)')
# #plt.savefig('u_0.pdf',bbox_inches='tight')

# uvel = mit.rdmds(dir0 + file1,-1)
# plt.figure()
# plt.contourf(1e-3*YC[:,np.int(si_x/2)],RC[:,0,0],uvel[:,:,np.int(si_x/2)],30,cmap=plt.cm.seismic)
# plt.colorbar(label=r'u (m\,s$^{-1}$)')
# plt.xlabel('y(km)')
# plt.ylabel('z(m)')
# #plt.savefig('u_1.pdf',bbox_inches='tight')


# for i in iters:
#   data = mit.rdmds(dir0 + file1,i)
#   psi = data[50,:,:]
# #  psi[0,0] = -0.15
# #  psi[0,1] =  0.15
#   plt.figure()
#   plt.contourf(XC*1e-3,YC*1e-3,psi,40,cmap=plt.cm.seismic)
#   plt.xlabel('km')
#   plt.ylabel('km')
#   plt.colorbar(label='m/s')
