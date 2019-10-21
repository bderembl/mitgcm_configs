#!/usr/bin/env python

import numpy as np
import scipy.io.netcdf as netcdf

binprec = '>f4'

def read_var(dir0, var, app):
  f1 = netcdf.netcdf_file(dir0 + var + '.' + app,'r',maskandscale=True)
  psi = f1.variables[var]
  psi = psi[:,::-1,:]
  f1.close()
  return psi

# ERA 40 data files (cf. https://apps.ecmwf.int/datasets/data/era40-daily/levtype=sfc/)
# you can download a sample here: 
# www.lmd.ens.fr/deremble/mitgcm_configs/cheapaml_era40/data/

dir0 = '../data/'
date = ['200001']

nfi = 0
# ================ OCEAN TEMPERATURES =========================================

var = 'sstk'
app = '.ash.nc'

sst = read_var(dir0, 'sstk', date[nfi] + '.ash.nc')

# in degC
sst = sst - 273.16
# land and sea ice
sst = np.where(sst>350,0.,sst)
sst = np.where(sst<0,0.,sst)

sst.astype(binprec).tofile('temp.box')

# ======================= GRID    ========================================== 
si_t,si_y,si_x = sst.shape
si_z = 1

dx = 360/si_x*np.ones(si_x) # 360 deg long
dy = 170/si_y*np.ones(si_y) # 170 deg lat (not 180 to avoid singularity at the pole)
dz = 100.*np.ones(si_z)     # arbitrary

dx.astype(binprec).tofile('dx.box')
dy.astype(binprec).tofile('dy.box')
dz.astype(binprec).tofile('dz.box')

# =======================  TOPO AND MASK     ========================================== 
h = -dz*np.ones((si_y,si_x))

# adjust land sea mask
h = np.where(sst<=0,0.,h)
mask = np.where(h<0,1.,0)

h.astype(binprec).tofile('topog.box')
mask.astype(binprec).tofile('rbcs_mask.box')


# ================ ATMOSPHERIC TEMPERATURES =========================================

t2m = read_var(dir0, 't2', date[nfi] + '.ash.nc')
t2m = t2m - 273.16

t2m.astype(binprec).tofile('tair.box')

# ================ ATMOSPHERIC HUMIDITY =========================================

d2m = read_var(dir0, 'd2', date[nfi] + '.ash.nc')
d2m = d2m - 273.15
# convert dew point to specific humidity (kg/kg) at 1000 mb
qair = 6.112*np.exp(17.67*d2m/(243.5+d2m))*0.622/1000

qair.astype(binprec).tofile('qair.box')

# =======================  WIND     ========================================== 

u10 = read_var(dir0, 'u10', date[nfi] + '.ash.nc')
u10.astype(binprec).tofile('windx.box')

v10 = read_var(dir0, 'v10', date[nfi] + '.ash.nc')
v10.astype(binprec).tofile('windy.box')

# =======================  SOLAR     ========================================== 

ssrd = read_var(dir0, 'ssrd', date[nfi] + '.fsh.nc')

# in W/m2
ssrd = ssrd/60./60./6.

ssrd.astype(binprec).tofile('solar.box')

# =======================  CLOUD FRACTION     ========================================== 

tcc = read_var(dir0, 'tcc', date[nfi] + '.ash.nc')
tcc.astype(binprec).tofile('cloud.box')

# ======================= BOUNDARY LAYER   ========================================== 

blh = read_var(dir0, 'blh', date[nfi] + '.fsh.nc')
blh.astype(binprec).tofile('h_blayer.box')

print("End time and forcing cycle: {}".format(si_t*6*3600))
