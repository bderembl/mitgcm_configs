#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
from scipy import interpolate

binprec = '>f4'
plt.ion()

flag_buildtopo = 1

#% ================ GRID =========================================                                  

si_x = 300
si_y = 150
  
#nb of layers
si_z = 1

# max depth
Lz = 1000.0


lat1 = 15
lat2 = 60

lon1 = -100
#lon2 = 28
lon2 = 0

R    = 6400e3
deg1 = 2*np.pi*R/360

lat = np.linspace(lat1,lat2,si_y)
lon = np.linspace(lon1,lon2,si_x)
lon_g,lat_g = np.meshgrid(lon,lat) 


xx = np.zeros((si_x*si_y,2))
xx[:,0] = lat_g.flat
xx[:,1] = lon_g.flat

dx = lon[1] - lon[0]
dy = lat[1] - lat[0]

dx1_m = dx*deg1*np.abs(np.cos(lat1*np.pi/180.))
dx2_m = dx*deg1*np.abs(np.cos(lat2*np.pi/180.))

print('dx1 = {0:.1f} m'.format(dx1_m))
print('dx2 = {0:.1f} m'.format(dx2_m))
print('dy = {0:.1f} m'.format(dy*deg1))

dx1 = dx*np.ones((si_x))
dy1 = dy*np.ones((si_y))
dz1 = Lz*np.ones(1)

dx1.astype(binprec).tofile( 'dx.box')
dy1.astype(binprec).tofile( 'dy.box')
dz1.astype(binprec).tofile( 'dz.box')

print('config: si_x = {0}, si_y = {1}, si_z = {2}'.format(si_x,si_y,si_z))

# ===== TOPOGRAPHY =======
if flag_buildtopo:
  
  file_topo = '/home/bderembl/work/data/topo/etopo1.nc'
  lat_name = 'y'
  lon_name = 'x'
  topo_name = 'z'
  
  ft = netcdf.netcdf_file(file_topo,'r')
  
  lat_t = ft.variables[lat_name][:].copy().squeeze()
  lon_t = ft.variables[lon_name][:].copy().squeeze()
  
  
  y1 = np.argmin(np.abs(lat_t-lat1))-10
  y2 = np.argmin(np.abs(lat_t-lat2))+10
  
  x1 = np.argmin(np.abs(lon_t-lon1))-10
  x2 = np.argmin(np.abs(lon_t-lon2))+10
  
  topo =  ft.variables[topo_name][y1:y2,x1:x2].copy().squeeze()
  
  lat_t = lat_t[y1:y2]
  lon_t = lon_t[x1:x2]
  
  lon_t2,lat_t2 = np.meshgrid(lon_t,lat_t) 
  
  func_topo = interpolate.interp2d(lon_t, lat_t, topo, kind='linear')
  
  topo_newgrid = func_topo(lon,lat)
  
  topo_newgrid = np.where(topo_newgrid > 0, 0, topo_newgrid)
  topo_newgrid = np.where(topo_newgrid < 0, -Lz, topo_newgrid)

  # Walls
  topo_newgrid[:,0] = 0.0
  topo_newgrid[0,:] = 0.0
  
  topo_newgrid.astype(binprec).tofile('htopo.box')
  

# ===== TOPOGRAPHY =======

windx = np.zeros((2,si_y,si_x))

# stess
windx[0,:,:] = 1e-5*(lat_g - lat1)/(lat2-lat1)
windx[1,:,:] = -windx[0,:,:]

windx.astype(binprec).tofile('windx.box')
