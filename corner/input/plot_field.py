#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf

plt.ion()

flag_mov = 0
flag_traj = 0

dir0 = '../run/'

file1 = 'diags.0000000000.t001.nc'
file2 = 'grid.t001.nc'

f1 = netcdf.netcdf_file(dir0 + file1)
f2 = netcdf.netcdf_file(dir0 + file2)


x = f2.variables['X'][:].copy()
y = f2.variables['Y'][:].copy()

xp1 = f2.variables['Xp1'][:].copy()
yp1 = f2.variables['Yp1'][:].copy()
T = f1.variables['T'][:].copy()


si_x = len(x)
si_y = len(y)
si_t = len(T)

h_mit = f2.variables['Depth'][:,:].copy()

vort = f1.variables['momVort3'][0,:,:].copy()

vmin = np.min(vort)
vmax = -vmin
vcont = np.linspace(vmin,vmax,20)


xunit = 1000.0 # 1:m -- 1000:km

posxy = np.zeros((2,si_t),dtype='int')

if flag_traj == 1:
  for nt in range(0,si_t):
    vort = f1.variables['momVort3'][nt,:,:].copy()
    posxy[0,nt],posxy[1,nt] = np.unravel_index(np.argmin(vort),vort.shape)
    

plt.figure()

if flag_mov == -1:
  nt = 0
  mytime = [49]
  vort = f1.variables['momVort3'][mytime[nt],:,:].copy()
  plt.contour(xp1[:si_x/2]/xunit,yp1/xunit,vort[:,:si_x/2],vcont,colors='k')
  plt.title('Day ' + str(mytime[nt]+1))
  plt.xlabel('x (km)')
  plt.ylabel('y (km)')
  myci = "CI: {:.1e}".format(vcont[1]-vcont[0])
  plt.text(x[120]/xunit,y[5]/xunit,myci)

  if flag_traj:
    plt.plot(xp1[posxy[1,:mytime[nt]]]/xunit,yp1[posxy[0,:mytime[nt]]]/xunit,'b')
    plt.plot(xp1[posxy[1,mytime[nt]:]]/xunit,yp1[posxy[0,mytime[nt]:]]/xunit,'b--')

elif flag_mov == 0:
  mytime = [0,9,19,29]

  for nt in range(0,len(mytime)):
    plt.subplot(2,2,nt+1, aspect='equal')
    vort = f1.variables['momVort3'][mytime[nt],:,:].copy()
    plt.contour(xp1/xunit,yp1/xunit,vort.squeeze(),vcont,colors='k')
    plt.contourf(x/xunit,y/xunit,h_mit,[-10,0],colors='0.5')
    plt.title('Day ' + str(mytime[nt]+1))
    if nt == 2 or nt == 3:
      plt.xlabel('x (km)')
    if nt == 0 or nt == 2:      
      plt.ylabel('y (km)')
    myci = "CI: {:.1e}".format(vcont[1]-vcont[0])
    plt.text(x[-170]/xunit,y[5]/xunit,myci)

  plt.savefig('corner_10mit.eps')

elif flag_mov == 1:

  vort = f1.variables['momVort3'][:,:,:].copy()

  vmin = np.min(vort)
  vmax = -vmin

  vcont = np.linspace(vmin,vmax,20)

  for nt in range(0,si_t):
    vort = f1.variables['momVort3'][nt,:,:].copy()
    vort = vort.squeeze()
    vort[0,0] = vmin
    vort[0,1] = vmax
    plt.contourf(xp1/xunit,yp1/xunit,vort,vcont,cmap = plt.cm.bwr)
    plt.contourf(x/xunit,y/xunit,h_mit,[-10,0],colors='0.5')
    ext = '0'
    if nt > 9:
      ext = ''
    plt.savefig('movie/ewall_'+ ext + str(nt) + 'mit.png') 
    plt.clf()

f1.close()
f2.close()
