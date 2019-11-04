#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit

dir0 = '/run/media/bderembl/izar/data/test_pg_hr/run06/'

file_u = 'U*'
file_v = 'V*'
file_w = 'W*'
file_t = 'T*'
file_s = 'S*'
file_d = 'oceDiag*'

flag_grid = 1
flag_comp_mean = 1 #1: compute and save mean, 0: read from file


# load grid
if flag_grid:
  XC    = mit.rdmds(dir0+'XC*')
  YC    = mit.rdmds(dir0+'YC*')
  XG    = mit.rdmds(dir0+'XG*')
  YG    = mit.rdmds(dir0+'YG*')
  DXC   = mit.rdmds(dir0+'DXC*')
  DYC   = mit.rdmds(dir0+'DYC*')
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

iters1 = mit.mds.scanforfiles(dir0 + file_u)
iters2 = mit.mds.scanforfiles(dir0 + file_d)

si_t1 = len(iters1)
si_t2 = len(iters2)

i = -1
uvel = mit.rdmds(dir0 + file_u,iters1[i],rec=0)
si_z,si_y,si_x = uvel.shape

temp = mit.rdmds(dir0 + file_t,iters1[i],rec=0)

if flag_comp_mean:
  gT_Forc_me = np.zeros((si_z,si_y,si_x))
  for it in iters2:
    gT_Forc_me = gT_Forc_me + mit.rdmds(dir0 + file_d,it,rec=0)

  gT_Forc_me = gT_Forc_me/si_t2
  gT_Forc_me.tofile('gT_Forc_me.box')
else:
  gT_Forc_me = np.fromfile('gT_Forc_me.box',dtype=float, count=-1, sep='').reshape(si_z,si_y,si_x)

sc = 1e-3
vcont = np.linspace(-2e-6,2e-6,50)

plt.figure()
plt.contourf(sc*XC,sc*YC,gT_Forc_me[1,:,:],vcont,cmap=plt.cm.bwr,extend='both')

