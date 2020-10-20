#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
from scipy.interpolate import interp1d
import scipy.special, scipy.interpolate, scipy.integrate

plt.ion()

binprec = '>f8'

flag_out = 0 #0: mitgcm restart, 1: qg, 2: HIM

#% ================== NEW GRID =====================================
rSphere = 6370.e3
deg2m = 2*np.pi*rSphere/360.0


if flag_out == 0: # mitgcm
#  si_x = 100
#  si_y = 100
  si_x = 500
  si_y = 500
  Lx = 500.0e3
  Ly = 500.0e3
elif flag_out == 1: # q-gcm
  si_x = 481
  si_y = 481
  # in m
  Lx = 500.0e3
  Ly = 500.0e3
elif flag_out == 2: # HIM
  si_x = 500
  si_y = 500
  # in m
  Lx = 500.0e3
  Ly = 500.0e3

si_z = 100;


si_x1 = si_x + 1
si_y1 = si_y + 1

Lz = 1000

dx = Lx/si_x;
dy = Ly/si_y;
dz = Lz/si_z

xx = Lx*(np.arange(0,si_x) + 0.5)/(1.0*si_x)
yy = Ly*(np.arange(0,si_y) + 0.5)/(1.0*si_y)

xx1 = Lx*(np.arange(0,si_x+1) )/(1.0*si_x)
yy1 = Ly*(np.arange(0,si_y+1) )/(1.0*si_y)


xg,yg = np.meshgrid(xx,yy) 
xu,yu = np.meshgrid(xx1[:-1],yy) 
xv,yv = np.meshgrid(xx,yy1[:-1]) 
xc,yc = np.meshgrid(xx1,yy1) 


dx1 = dx*np.ones((si_x))
dy1 = dy*np.ones((si_y))
dz1 = dz*np.ones((si_z))


zz = dz1.cumsum()
zc = np.reshape(zz - 0.5*dz1,(si_z,1,1))
zf = np.zeros(si_z+1)
zf[1:] = zz
zf = np.reshape(zf,(si_z+1,1,1))

dx1.astype(binprec).tofile('dx.box')
dy1.astype(binprec).tofile('dy.box')
dz1.astype(binprec).tofile('dz.box')



#% ============== background density profile ===================

temp_i  = 5.0*(1 - zc/Lz)
temp_if = 5.0*(1 - zf/Lz)

#%==================== SST - LAND ===================================

landh  = np.zeros((si_y,si_x));
theta  = np.zeros((si_z,si_y,si_x));
uvel   = np.zeros((si_z,si_y,si_x));
vvel   = np.zeros((si_z,si_y,si_x));


H = dz1.cumsum()[-1]
landh = -H + landh

landh[:,0] = 0.0
# wall at northern and southern boundaries
landh[0,:] = 0.0
landh[-1,:] = 0.0
landh.astype(binprec).tofile('topog.box')


# # physical constants
#rho_const = 999.8
rho_const = 1000.0
alphaK = 2.0e-4
#g0 = 9.8
g0 = 10.0
f0 = 1.0e-4

# 1: lamb-oseen, 2: rankine
flag_vortex_type = 2


#gaussian eddy
x_c = 75.0e3
y_c = 200.0e3
DeltaT = 5.0
if flag_vortex_type == 1:
  # lamb-oseen
  R0 = 30e3
  #Gamma (lamb-oseen): cyclone: -1 anticyclone: +1
  Gamma = +1.0e5
elif flag_vortex_type == 2:
  # rankine
  R0 = 40e3
  #vmax (rankine): cyclone: +1 anticyclone: -1
  vmax = 1.0

z0 = 500.0 # m



# ERF vertical profile
FZ = 1-scipy.special.erf(zc/z0)
# vertical derivative
FpZ  = -1/z0*2/np.sqrt(np.pi)*np.exp(-zc**2/z0**2)
FpZf = -1/z0*2/np.sqrt(np.pi)*np.exp(-zf**2/z0**2)


FZ   = FZ.reshape(si_z,1,1)
FpZ  = FpZ.reshape(si_z,1,1)
FpZf = FpZf.reshape(si_z+1,1,1)


# grid at U,V,T points
rad_gg = np.sqrt((xg-x_c)**2 + (yg-y_c)**2)
rad_cc = np.sqrt((xc-x_c)**2 + (yc-y_c)**2)
rad_gu = np.sqrt((xu-x_c)**2 + (yu-y_c)**2)
rad_gv = np.sqrt((xv-x_c)**2 + (yv-y_c)**2)

theta_gg = np.arctan2(yg-y_c,xg-x_c)
theta_cc = np.arctan2(yc-y_c,xc-x_c)
theta_gu = np.arctan2(yu-y_c,xu-x_c)
theta_gv = np.arctan2(yv-y_c,xv-x_c)

# image effect
rad_gg_i = np.sqrt((xg+x_c)**2 + (yg-y_c)**2)
rad_cc_i = np.sqrt((xc+x_c)**2 + (yc-y_c)**2)
rad_gu_i = np.sqrt((xu+x_c)**2 + (yu-y_c)**2)
rad_gv_i = np.sqrt((xv+x_c)**2 + (yv-y_c)**2)

theta_gg_i = np.arctan2(yg-y_c,xg+x_c)
theta_cc_i = np.arctan2(yc-y_c,xc+x_c)
theta_gu_i = np.arctan2(yu-y_c,xu+x_c)
theta_gv_i = np.arctan2(yv-y_c,xv+x_c)


# create velocity profile
# lamb-oseen vortex tangent velocity profile
def vel_lamb(rr):
  v = Gamma/(2*np.pi*rr)*(1-np.exp(-rr**2/R0**2))
  v = np.where(rr == 0, 0.0,v)
  return v

# hyperpolic vortex (~ rankine)
def vel_rankine(rr):
  v = -vmax*np.tanh(rr/R0)/(np.cosh(rr/R0))**2
  v = np.where(rr == 0, 0.0,v)
  return v

if flag_vortex_type == 1:
  u_out = vel_lamb(rad_gu)*np.sin(-theta_gu)
  v_out = vel_lamb(rad_gv)*np.cos(theta_gv)

  u_out_i = -vel_lamb(rad_gu_i)*np.sin(-theta_gu_i)
  v_out_i = -vel_lamb(rad_gv_i)*np.cos( theta_gv_i)
  # 3D vorticity profile
  vort =  Gamma/(np.pi*R0**2)*(np.exp(-rad_cc**2/R0**2))
elif flag_vortex_type == 2:
  u_out = vel_rankine(rad_gu)*np.sin(-theta_gu)
  v_out = vel_rankine(rad_gv)*np.cos(theta_gv)

  u_out_i = -vel_rankine(rad_gu_i)*np.sin(-theta_gu_i)
  v_out_i = -vel_rankine(rad_gv_i)*np.cos( theta_gv_i)


# real + image
u_out = u_out + u_out_i
v_out = v_out + v_out_i

# 3D velocity field
uvel = FZ*np.tile(u_out,[si_z,1,1])
vvel = FZ*np.tile(v_out,[si_z,1,1])


# compute pressure field
def integrand(rr):
  if flag_vortex_type == 1:
    v = vel_lamb(rr)
  elif flag_vortex_type == 2:
    v = vel_rankine(rr)

# neglect cyclostrophic for now to keep it linear
  res = v**2/rr  + f0*v
#  res = f0*v
  res = np.where(rr == 0, 0.0,res)
  return res

def comp_p(x):
  
  if x ==0:
    xx = 1e-12
  else:
    xx = 1.0*x

  a,b = scipy.integrate.quad(integrand,0,xx)
  return a

Lmax = np.max([Lx,Ly])
si_max = np.max([si_x,si_y])
rr = np.linspace(0.0,2*Lmax, 10*si_max)
p1 = [ comp_p(x) for x in rr.flatten() ]
fint = scipy.interpolate.interp1d(rr, p1)

p_out   =  rho_const*fint(rad_gg)
p_out_i = -rho_const*fint(rad_gg_i)

p_out = p_out + p_out_i
# remove const at infinity
p_out = p_out - p_out[0,0]


eta = p_out/rho_const/g0
pres = FZ*np.tile(p_out,[si_z,1,1])

# minus sign
#dpdz = np.diff(pres,1,0)/dz2
dpdz = FpZ*np.tile(p_out,[si_z,1,1])
rhop = dpdz/g0

# convert to temperature
theta_a = -rhop/(rho_const*alphaK) 
theta = theta_a + temp_i

# face (with surface)
rhopf = FpZf*np.tile(p_out,[si_z+1,1,1])/g0

# convert to temperature
theta_af = -rhopf/(rho_const*alphaK) 
thetaf = theta_af + temp_if

uvel.astype(binprec).tofile('uinit.box')
vvel.astype(binprec).tofile('vinit.box')
theta.astype(binprec).tofile('tinit.box')

eta.astype(binprec).tofile('einit.box')

temp_i.astype(binprec).tofile('tref.box')



# print infos
print("barotropic Rd: {} km".format(np.sqrt(H*g0)/f0/1000.0))

if flag_out == 1:
  # prepare a QG restart file
  # Store 
  fileout = 'lastday_r.nc'
  f = netcdf.netcdf_file('/home/bderembl/work/bqg/q-gcm/examples/eddy_wall/' + fileout,'w')
  
  si_zq = 10
  nbc = 2*(si_y-3) + 2*(si_x-3) + 2
  
  f.createDimension('ypo',si_y)
  f.createDimension('xpo',si_x)
  f.createDimension('zo' ,si_zq)
  f.createDimension('yto',si_y-1)
  f.createDimension('xto',si_x-1)
  f.createDimension('xbc',nbc)
  f.createDimension('time',1)
  
  ypo = f.createVariable('ypo', 'd', ('ypo',))
  xpo = f.createVariable('xpo', 'd', ('xpo',))
  zo  = f.createVariable('zo' , 'd', ('zo' ,))
  yto = f.createVariable('yto', 'd', ('yto',))
  xto = f.createVariable('xto', 'd', ('xto',))
  xbc = f.createVariable('xbc', 'd', ('xbc',))
  time = f.createVariable('time', 'd', ('time',))
  
  
  po  = f.createVariable('po' , 'd', ('zo','ypo','xpo',))
  pom  = f.createVariable('pom' , 'd', ('zo','ypo','xpo',))
  
  sst   = f.createVariable('sst' , 'd', ('yto','xto',))
  sstm  = f.createVariable('sstm' , 'd', ('yto','xto',))

  vw  = f.createVariable('vw' , 'd', ('zo','xbc',))
  Mw  = f.createVariable('Mw' , 'd', ('zo','xbc',))
  
  

  ypo[:] = np.arange(si_y)
  xpo[:] = np.arange(si_x)
  zo[:]  = np.arange(si_zq)
  yto[:] = np.arange(si_y-1)
  xto[:] = np.arange(si_x-1)
  xbc[:] = np.arange(nbc)
  time[:] = 0.0
  
  po [:,:,:] = pres[::10,:,:]/rho_const
  pom [:,:,:] = pres[::10,:,:]/rho_const
  sst [:,:]  = 0.0*theta[0,1:,1:]
  sstm [:,:] = 0.0*theta[0,1:,1:]
  
  vw[:,:] = np.zeros((si_zq,nbc))
  Mw[:,:] = np.zeros((si_zq,nbc))

  f.close()

if flag_out == 2:
  # prepare a GOLD restart
  si_z2 = 10
  gfs  = 0.98
  gp = 10*1e3*2e-4*DeltaT/si_z2/1e3

  # background layer thickness
  h0 = Lz/si_z2
  
  # hi: interface
  hi = np.zeros((si_z2+1,si_y,si_x))

  idzskip = np.int(si_z/si_z2)
  pres2 = pres[::idzskip,:,:]
  etag = pres[0,:,:]/rho_const/gfs
  hi[0,:,:] = etag
  
  # method 1
  
  for nz in range(1,si_z2):
    hi[nz,:,:] = (pres2[nz,:,:] - pres2[nz-1,:,:])/rho_const/gp  
  ht = -np.diff(hi,1,0) + h0

  # derive geostrophic velocities
  u = -1.0/f0/rho_const*np.diff(pres2,1,1)/dy
  v =  1.0/f0/rho_const*np.diff(pres2,1,2)/dx
  
  u1 = np.zeros((si_z2,si_y,si_x1))
  v1 = np.zeros((si_z2,si_y1,si_x))
  
  u1[:,:-1,:-1] = u
  v1[:,:-1,:-1] = v

  # # method 2: with outcrop (si_z2 is unknowm)
  rhofull = -thetaf*rho_const*alphaK
  # set bottom rho to zero for interpolation
  rhofull[-1,:,:] = 1e-8
  
  rhomin = np.min(rhofull)
  rhomax = np.max(rhofull)

  drho = gp/g0*rho_const
  if DeltaT > 0:
    rhot = np.arange(0.0,-rhomin,drho)
    rhot = -rhot[::-1]
    # set min rho (merge upper 2 layers)
    rhofull = np.where(rhofull<=rhot[0],rhot[0]+1e-5,rhofull)
  else:
    # do something
    print("write cold anomaly case")


    
  si_z2 = len(rhot) - 1
  hi = np.zeros((si_z2+1,si_y,si_x))
  for nx in range(0,si_x):
    for ny in range(0,si_y):

      fint = scipy.interpolate.interp1d(rhofull[:,ny,nx],zf.squeeze(),bounds_error=False,fill_value=-etag[ny,nx])
      hi[:,ny,nx] = fint(rhot)

      
  #merge upper two layers
  hi = -hi
  
  ht = -np.diff(hi,1,0)

  # add 0.1 m to each layer
  ht = ht + 0.1
  
  pres_lay = 0.0*ht
  
  pres_lay[0,:,:] = hi[0,:,:]*gfs
  for k in range(1,si_z2):
    pres_lay[k,:,:] = gp*hi[k,:,:]
    
  pres_lay = np.cumsum(pres_lay,0)

  # derive geostrophic velocities
  u = -1.0/f0*np.diff(pres_lay,1,1)/dy
  v =  1.0/f0*np.diff(pres_lay,1,2)/dx
  
  u1 = np.zeros((si_z2,si_y,si_x1))
  v1 = np.zeros((si_z2,si_y1,si_x))
  
  u1[:,:-1,:-1] = u
  v1[:,:-1,:-1] = v

  # get velocities via interpolation
  zmid = -0.5*(hi[:-1,:,:] + hi[1:,:,:])
  for nx in range(0,si_x):
    for ny in range(0,si_y):

      fuint = scipy.interpolate.interp1d(zc.squeeze(),uvel[:,ny,nx].squeeze(),fill_value="extrapolate")
      fvint = scipy.interpolate.interp1d(zc.squeeze(),vvel[:,ny,nx].squeeze(),fill_value="extrapolate")
      u1[:,ny,nx] = fuint(zmid[:,ny,nx])
      v1[:,ny,nx] = fvint(zmid[:,ny,nx])

      
  
  # -- write netcdf --
  
  file1 = 'myrestart.nc'
  dir_him = '/home/bderembl/work/HIM/myrun/eddy-wall/saves/'
  
  f1 = netcdf.netcdf_file(dir_him + file1,'w')
  
  f1.createDimension('Time',1)
  f1.createDimension('Layer',si_z2)
  f1.createDimension('lath',si_y)
  f1.createDimension('lonh',si_x)
  f1.createDimension('latq',si_y1)
  f1.createDimension('lonq',si_x1)
  
  time_o  = f1.createVariable('Time','f',('Time',))
  layer_o = f1.createVariable('Time','f',('Layer',))
  lath_o  = f1.createVariable('lath','f',('lath',))
  lonh_o  = f1.createVariable('lonh','f',('lonh',))
  latq_o  = f1.createVariable('latq','f',('latq',))
  lonq_o  = f1.createVariable('lonq','f',('lonq',))
  
  u_o  = f1.createVariable('u','f',('Time','Layer','lath','lonq',))
  v_o  = f1.createVariable('v','f',('Time','Layer','latq','lonh',))
  h_o  = f1.createVariable('h','f',('Time','Layer','lath','lonh',))
  tr_o = f1.createVariable('tr0','f',('Time','Layer','lath','lonh',))
  
  time0 = 0.0
  dep  = np.arange(si_z2)
  lath = np.arange(si_y)
  lonh = np.arange(si_x)
  latq = np.arange(si_y1)
  lonq = np.arange(si_x1)
  
  time_o[:]  = time0
  layer_o[:] = dep
  lath_o[:]  = lath
  lonh_o[:]  = lonh
  latq_o[:]  = latq
  lonq_o[:]  = lonq
  
  u_o [:,:,:,:] = u1
  v_o [:,:,:,:] = v1
  h_o [:,:,:,:] = ht
  tr_o[:,:,:,:] = ht
  
  
  f1.close()
  
  
  # -- write netcdf (D) --
  file2 = 'MESO_topog_2.cdf'
  
  f2 = netcdf.netcdf_file(dir_him + file2,'w')
  
  f2.createDimension('lath',si_y)
  f2.createDimension('lonh',si_x)
  
  lath_o  = f2.createVariable('lath','f',('lath',))
  lonh_o  = f2.createVariable('lonh','f',('lonh',))
  
  D_o  = f2.createVariable('D','f',('lath','lonh',))
  
  lath_o[:]  = lath
  lonh_o[:]  = lonh
  
  D_o[:,:]  = H
  
  f2.close()



  lenlon = Lx/111317.0
#  print lenlon
