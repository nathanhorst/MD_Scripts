import sys
import os
import numpy as np

import add_poly as ap

#################
## File to read
#################

inputfile='Au201.xyz'
outputfile='Au201_halfPEO.xyz'

#############################################
# Define coordinates of points to build from
#############################################\

x1=0.0
y1=0.0
z1=2.0

R1=np.sqrt(x1**2+y1**2+z1**2)

shift1=0.5

###################################
# Length of polymer chain to build
###################################

L=10

##############
# Run Builder
##############

ap.poly_without_points(inputfile,outputfile,R1,L,shift1)
