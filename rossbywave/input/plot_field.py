#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
plt.ion()

flag_mov = 1

dir0 = '../run/mnc_test_0002/'

file1 = 'state.0000000000.t001.nc'
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

ssh = f1.variables['Eta'][-1,:,:].copy()

vmin = np.min(ssh)
vmax = -vmin
vcont = np.linspace(vmin,vmax,20)


xunit = 1000.0 # 1:m -- 1000:km


plt.figure()

if flag_mov == -1:
  nt = 0
  mytime = [49]
  psi = f1.variables['Eta'][mytime[nt],:,:].copy()
  ax = plt.subplot(projection=ccrs.PlateCarree());
  ax.contourf(x, y, psi, vcont,transform=ccrs.PlateCarree(),cmap='bwr')
  ax.coastlines()
  gl = ax.gridlines(draw_labels=True, alpha = 0.5, linestyle='--');
  gl.xlabels_top = False
  gl.ylabels_right = False
  gl.xformatter = LONGITUDE_FORMATTER
  gl.yformatter = LATITUDE_FORMATTER

  plt.title('Day ' + str(mytime[nt]+1))

elif flag_mov == 0:
  mytime = [0,9,19,29]

  for nt in range(0,len(mytime)):
    psi = f1.variables['Eta'][mytime[nt],:,:].copy()

    ax = plt.subplot(2,2,nt+1, aspect='equal',projection=ccrs.PlateCarree());
    ax.contourf(x, y, psi, vcont,transform=ccrs.PlateCarree(),cmap='bwr')
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True, alpha = 0.5, linestyle='--');
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    plt.title('Day ' + str(mytime[nt]+1))

  plt.savefig('rossby.eps')

elif flag_mov == 1:

  psi = f1.variables['Eta'][:,:,:].copy()

  vmin = np.min(psi)
  vmax = -vmin

  vcont = np.linspace(vmin,vmax,20)

  for nt in range(0,si_t):
    psi = f1.variables['Eta'][nt,:,:].copy().squeeze()

    psi[0,0] = vmin
    psi[0,1] = vmax

    ax = plt.subplot(projection=ccrs.PlateCarree());
    ax.contourf(x, y, psi, vcont,transform=ccrs.PlateCarree(),cmap='bwr')
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True, alpha = 0.5, linestyle='--');
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    plt.title('Day ' + str(T[nt]/86400))


    ext = '00'
    if nt > 9:
      ext = '0'
    if nt > 99:
      ext = ''
    plt.savefig('movie/ssh_'+ ext + str(nt) + 'mit.png') 
    plt.clf()

f1.close()
f2.close()
