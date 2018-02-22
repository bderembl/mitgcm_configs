#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import scipy.io.netcdf as netcdf
from scipy.special import erf

plt.ion()

dir1 = '../run/mnc_test_0001/'
file1 = 'state.0000000000.t001.nc'


f1 = netcdf.netcdf_file(dir1 + file1,'r')

z = f1.variables['Z'][:].copy()
time = f1.variables['T'][:].copy()
u = f1.variables['U'][:,:,0,0].copy()


plt.figure()

it = -1
plt.plot(u[it,:],z,'r',label='model')
plt.plot(erf((z-z[-1])/(2*np.sqrt(time[it]))),z,'k--',label='similarity solution')
plt.xlabel('U')
plt.ylabel('z')
plt.legend()
plt.savefig('similarity_sol.png')
