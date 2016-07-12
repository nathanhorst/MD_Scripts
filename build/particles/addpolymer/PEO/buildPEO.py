# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 14:43:52 2016

@author: nathan
"""
import numpy as np

import build_PEO_np as bpnp

#############################
## Number of repeating units in polymer graft
#############################

n=20

#################
## File to output grafted nanoparticle
#################

outputfile1='Au201_PEOhexchains.xyz' 

#################
## xyz File with existing nanoparticle
#################

npfile1='Au201shell.xyz'

#################
## File to store polymer before grafting
#################

poly_file='PEO_chain_n'+str(n)+'.xyz'

#################
## build and save polymer
#################

bpnp.build_PEO_chain_noS(poly_file,n)

#################
##set grafting vectors
#################

x1=1
y1=1
z1=-1

R1=np.sqrt(x1**2+y1**2+z1**2)

###########

x2=-1
y2=0.5
z2=-1.5

R2=np.sqrt(x2**2+y2**2+z2**2)

#################
## graft chain and save nanoparticle
#################

bpnp.graft_PEO_hex(outputfile1,poly_file,npfile1,R1,R2)

#################
## File to output grafted nanoparticle
#################

outputfile2='Au201PEO.xyz' 

#################
## xyz File with existing nanoparticle
#################

npfile2='Au201_PEOhexchains.xyz'

#################
##set grafting vectors
#################

x1=0
y1=0
z1=2

R1=np.sqrt(x1**2+y1**2+z1**2)

###########
##s-s is the sulfur-sulfur distance on the square face in whatever units are being used
#############
s_s=4.1/3.96

#################
## graft chain and save nanoparticle
#################

bpnp.graft_square_faces(outputfile2,poly_file,npfile2,R1,s_s)

