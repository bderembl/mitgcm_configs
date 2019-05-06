#!/usr/bin/env python

import MITgcmutils as mit
import f90nml

dir0 = './'

fileA = 'pickup.ckptA.meta';
fileB = 'pickup.ckptB.meta';

filedata = 'data';

data = f90nml.read(dir0 + filedata)

ckptA = mit.mds.parsemeta(dir0 + fileA)
ckptB = mit.mds.parsemeta(dir0 + fileB)

if ckptA['timeStepNumber'] > ckptB['timeStepNumber']:
  ckpt = ckptA
  pickupsuff = 'ckptA'
else:
  ckpt = ckptB
  pickupsuff = 'ckptB'


patch_nml = {'parm03' : {'niter0' : ckpt['timeStepNumber'], 'pickupsuff' : pickupsuff}}

f90nml.patch('data', patch_nml,'newdata')
