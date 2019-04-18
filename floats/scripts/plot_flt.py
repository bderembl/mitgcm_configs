#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import time,datetime
import MITgcmutils
import glob,os,re

plt.ion()

dir1 = '../run/'

filef = 'float_trajectories*'

allfilesf = sorted(glob.glob(dir1 + filef));
si_f = len(allfilesf)
for nf in range(0,si_f):
  allfilesf[nf] = allfilesf[nf][:-5]

allfilesf = allfilesf[::2]
si_f = int(si_f/2)

# ==================== load floats ===============================

nf = 0
mit_fldata = MITgcmutils.rdmds(allfilesf[nf]).squeeze()
npy = int(allfilesf[nf][-3:])-1
npx = int(allfilesf[nf][-7:-4])-1

if si_f == 1:
  n1,n2 = mit_fldata.shape
  mit_fldata = mit_fldata.reshape(n1,n2,1,1)

dt = mit_fldata[0,1][npy,npx]

si_t  = mit_fldata[0,4][npy,npx].astype(int)
si_fl = mit_fldata[0,5][npy,npx].astype(int)

lat_fl = np.zeros((si_t,si_fl))
lon_fl = np.zeros((si_t,si_fl))


for nf in range(0,si_f):
  mit_fldata = MITgcmutils.rdmds(allfilesf[nf]).squeeze()

  npy = int(allfilesf[nf][-3:]) - 1
  npx = int(allfilesf[nf][-7:-4]) - 1

  if si_f == 1:
    n1,n2 = mit_fldata.shape
    mit_fldata = mit_fldata.reshape(n1,n2,1,1)

  n0 = 1
  
  flt_no = mit_fldata[1:,0][:,npy,npx].astype(int) - 1
  flt_it = (mit_fldata[1:,1][:,npy,npx]/dt).astype(int) - 1
  
  lon_flm = mit_fldata[1:,2][:,npy,npx]
  lat_flm = mit_fldata[1:,3][:,npy,npx]
  
  # lon_flm = lon_flm.reshape((si_t,si_fl))
  # lat_flm = lat_flm.reshape((si_t,si_fl))
  
  
  lon_fl[flt_it,flt_no] = lon_flm
  lat_fl[flt_it,flt_no] = lat_flm
  
nf = 100
ntmax = -1
for nf in range(0,si_fl-1):
  plt.plot(lon_fl[:,nf]*1e-3,lat_fl[:,nf]*1e-3,'r-',linewidth=1)
