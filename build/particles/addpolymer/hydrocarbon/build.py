# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 14:43:52 2016

@author: nathan
"""
import numpy as np

import build_HC_np as bhcnp

#############################
## Number of repeating units in polymer graft
#############################

n=11

#################
## File to output grafted nanoparticle
#################

outputfile1='Au201_HChexchains.xyz' 

#################
## xyz File with existing nanoparticle
#################

npfile1='Au201shell.xyz'

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

x1=3.96
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

bhcnp.graft_HC_hex(outputfile1,poly_file,npfile1,R1,R2)

#################
## File to output grafted nanoparticle
#################

outputfile2='Au201HC.xyz' 

#################
## xyz File with existing nanoparticle
#################

npfile2='Au201_HChexchains.xyz'

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

bhcnp.graft_square_faces(outputfile2,poly_file,npfile2,R1,s_s)

