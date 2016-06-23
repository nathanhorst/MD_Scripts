# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 14:43:52 2016

@author: nathan
"""
import os
import sys

import numpy as np
from numpy import linalg as LA

import build_HC_np as bhcnp
import build_HC_chain as bhcc

#############################
## Number of repeating units in polymer graft
#############################

n=11

#################
## File to output grafted nanoparticle
#################

outputfile='Au201_HCchains.xyz' 

#################
## xyz File with existing nanoparticle
#################

npfile='Au201_a3p96.xyz'

#################
## File to store polymer before grafting
#################

poly_file='HC_chain_n'+str(n)+'.xyz'

#################
## build and save polymer
#################

bhcnp.build_HC_chain_noS(poly_file,n)

#################
##set grafting vectors
#################

x1=-3.96
y1=3.96
z1=-3.96

R1=np.sqrt(x1**2+y1**2+z1**2)

###########

x2=-3.96
y2=1.98
z2=-5.94

R2=np.sqrt(x2**2+y2**2+z2**2)

#################
## graft chain and save nanoparticle
#################

bhcnp.graft_HC_chain(outputfile,poly_file,npfile,R1,R2)
