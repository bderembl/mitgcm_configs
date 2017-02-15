#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit

plt.ion()

dir0 = '/tank/groups/climode/bderembl/MITgcm/myrun/test_pg/run/'

dir1 = dir0 + 'mnc_test_0001/'
dir2 = dir0 + 'mnc_test*/'

file1 = 'state.0000000000.t*'


f1 = mit.mnc_files(dir2 + file1)



T = f1.variables['T'][:]
si_t = len(T)

nt = si_t-1
u = f1.variables['U'][nt,0,:,:]
temp = f1.variables['Temp'][nt,0,:,:]
eta = f1.variables['Eta'][nt,:,:]
