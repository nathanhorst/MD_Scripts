import sys
import os
import numpy as np

import add_poly as ap

#################
## File to read
#################

inputfile='Au201_halfPEO.xyz'
outputfile='Au201PEO.xyz'

#############################################
# Define coordinates of points to build from
#############################################

x1=1
y1=1
z1=1

R1=np.sqrt(x1**2+y1**2+z1**2)

###########

x2=1
y2=1.5
z2=0.5

R2=np.sqrt(x2**2+y2**2+z2**2)

###########
"""
x3=2
y3=3
z3=1

R3=np.sqrt(x3**2+y3**2+z3**2)
"""
###################################
# Length of polymer chain to build
###################################

L=10

##############
# Run Builder
##############

ap.poly_from_points(inputfile,outputfile,R1,R2,L)
