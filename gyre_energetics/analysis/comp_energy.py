#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import MITgcmutils as mit
import matplotlib

plt.ion()

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

dir0 = '../run/mnc_test_0001/'

files = 'state.0000000000.t001.nc'
fileg = 'grid.t001.nc'

# constants
rho0 = 1000
alphat = 2e-4
betas = 7.4e-4
gg = 10

# load grid
gridm = mit.rdmnc(dir0 + fileg)
XC    = gridm['XC'   ]
YC    = gridm['YC'   ]
XG    = gridm['XG'   ]
YG    = gridm['YG'   ]
DXC   = gridm['dxC'  ]
DYC   = gridm['dyC'  ]
hFacC = gridm['HFacC']
hFacS = gridm['HFacS']
hFacW = gridm['HFacW']
RAS   = gridm['rAs'  ]
RAW   = gridm['rAw'  ]
RAC   = gridm['rA'  ]
RAZ   = gridm['rAz'  ]
RC    = gridm['RC'   ]
RF    = gridm['RF'   ]
DRC   = gridm['drC'  ]
DRF   = gridm['drF'  ]
Depth = gridm['Depth']

# get dimensions and reshape
si_y, si_x = XC.shape
si_z, = RC.shape

si_x1 = si_x + 1
si_y1 = si_y + 1

XC  = np.reshape(XC,  (1, si_y,  si_x ))
YC  = np.reshape(YC,  (1, si_y,  si_x ))
XG  = np.reshape(XG,  (1, si_y1, si_x1))
YG  = np.reshape(YG,  (1, si_y1, si_x1))
RAC = np.reshape(RAC, (1, si_y,  si_x ))
DRF = np.reshape(DRF, (si_z, 1,  1 ))
RC  = np.reshape(RC,  (si_z, 1,  1 ))

# get 
iters = mit.rdmnc(dir0 + files, 'iter')['iter']

si_t = len(iters)

vel = mit.rdmnc(dir0 + files, ['U', 'V', 'T'], iters = [iters[-1]])
U = vel['U'].squeeze()
V = vel['V'].squeeze()
T = vel['T'].squeeze()

KE = 0.25*(U[:,:,:-1]**2 + U[:,:,1:]**2 + V[:,1:,:]**2 + V[:,:-1,:]**2)

ke_int = rho0*np.sum(KE*RAC*DRF)
#in PetaJoule
ke_int2 = ke_int*1e-15

b = gg*alphat*T

pe_int = rho0*np.sum(-b*RC*RAC*DRF)
#in PetaJoule
pe_int2 = pe_int*1e-15


print ('Total KE = {0:.2f} PJ'.format(ke_int2))
print ('Total PE = {0:.2f} PJ'.format(pe_int2))
