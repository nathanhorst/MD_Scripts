# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 15:52:26 2016

@author: waltmann
"""

import os
import sys

import numpy as np
from numpy import linalg as LA

import build_HC_np as bhcnp
import build_HC_chain as bhcc
import add_trans_poly as atp

#############################
## Number of repeating units in polymer graft
#############################

n=11

#################
## File to output grafted nanoparticle
#################

outputfile='complete.xyz' 

#################
## xyz File with existing nanoparticle
#################

npfile='Au201_HCchains.xyz'

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

x1=0
y1=0
z1=7.92

R1=np.sqrt(x1**2+y1**2+z1**2)

###########
##s-s is the sulfur-sulfur distance on the square face in whatever units are being used
#############
s_s=4.1

#################
## graft chain and save nanoparticle
#################


atp.graft_square_faces(outputfile,poly_file,npfile,R1,s_s)
