#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit

plt.ion()

#dir0 = '/home/bderembl/work/MITgcm/myrun/test_kw_energetics/run03/'
dir0 = '/media/bderembl/workd/MITgcm/myrun/test_kw_energetics/run08/'

dir1 = dir0 + 'mnc*/'
dir2 = dir0 + 'mnc_test_0001/'

file0 = 'grid.t*'
file1 = 'state.*'
file2 = 'oceDiag.*'

flag_mov = 1
flag_surf = 1

alphat = 2e-4
go = 9.81

# grid
f0 = mit.mnc_files(dir1 + file0)

RC    = f0.variables['RC'][:]
DRC   = f0.variables['drC'][:]
DRF   = f0.variables['drF'][:]
RF    = f0.variables['RF'][:]

XC    = f0.variables['XC'][:,:]
YC    = f0.variables['YC'][:,:]

si_y,si_x = XC.shape
si_z      = RC.size

dx = XC[1,1] - XC[0,0]
dy = YC[1,1] - YC[0,0]
dz = RC[1] - RC[0]

dv = np.abs(dx*dy*dz)

f1 = mit.mnc_files(dir1 + file1)
T = f1.variables['T'][:]
si_t = len(T)-1

f2 = mit.mnc_files(dir1 + file2)

nt = 0
if flag_mov:
  nt1 = 0
  nt2 = si_t
else:
  nt1 = nt
  nt2 = nt+1


# vort = f2.variables['momVort3'][:,0,:si_y,:si_x]
# vmin = 0.8*np.min(vort)
# vmax = -vmin
# vcont = np.linspace(vmin,vmax,50)

vmin = -5 #1.e-7
vmax = 5

vcont = 1e-5*np.linspace(vmin,vmax,51)
vmin = -5e-5
vmax = 5e-5

# section
nx_c = 75

plt.figure()

if flag_surf == 1:
  for nt in range(nt1,nt2):
    
    vort = f2.variables['momVort3'][nt,0,:si_y,:si_x]
      
    vort2 = np.where(vort>vmax,vmax,vort)
    vort2 = np.where(vort2<vmin,vmin,vort2)
  
    plt.contourf(XC*1e-3,YC*1e-3,vort2[:,:],vcont,cmap = plt.cm.bwr)
    #plt.plot([1e-3*XC[0,nx_c], 1e-3*XC[0,nx_c]],[1e-3*YC[0,nx_c], 1e-3*YC[-1,nx_c]],'k--')
    cb = plt.colorbar(label=r'$\omega$ (s$^{-1}$)')
  
    cb.formatter.set_powerlimits((0, 0))
    cb.update_ticks()
  
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
  
  
    
    app = '0'
    if nt > 9:
      app = ''
      
    if flag_mov:
      plt.savefig('figures/rv'+app+str(nt)+'.png',bbox_inches='tight')
      plt.clf()
    
    #plt.savefig('vort_surf.png',bbox_inches='tight')
  
  if flag_mov == 0:
    plt.figure()
  
else:
  # section at the wall
  #nx_c = 75
  nx_c = 1
  
  #temp = f1.variables['Temp'][:,:si_z,:si_y,nx_c]
  #v = f1.variables['V'][:,:si_z,:si_y,nx_c]
  
  vmin = 0.0#np.min(temp)
  vmax = 5.0#np.max(temp)
  vcont = np.linspace(vmin,vmax,11)
  
  vmin = 0.0#np.min(temp)
  vmax = 10.0#np.max(temp)
  vcont2 = np.linspace(vmin,vmax,21)
  
  vmin = -1.1
  vmax = 1.1
  vcont3 = np.linspace(vmin,vmax,51)
  
  for nt in range(nt1,nt2):
    
    
    temp = f1.variables['Temp'][nt,:si_z,:si_y,nx_c]
    v = f1.variables['V'][nt,:si_z,:si_y,nx_c]
  
  #  temp2 = np.where(temp>vmax,vmax,temp)
  #  temp2 = np.where(temp2<vmin,vmin,temp2) 
  
  #  plt.contourf(YC[:,0]*1e-3,RC,temp[:,:],vcont,cmap = plt.cm.viridis)
    plt.contourf(YC[:,0]*1e-3,RC,v[:,:],vcont3,cmap = plt.cm.seismic)
    plt.colorbar(ticks=[-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0])
    plt.contour(YC[:,0]*1e-3,RC,temp[:,:],vcont2,colors='k')
    #plt.savefig('temp_vsec.png',bbox_inches='tight')
  
    plt.xlabel('y (km)')
    plt.ylabel('z (m)')
  
    app = '0'
    if nt > 9:
      app = ''
      
    if flag_mov:
      plt.savefig('figures/twall_mit'+app+str(nt)+'.png',bbox_inches='tight')
      plt.clf()
    


