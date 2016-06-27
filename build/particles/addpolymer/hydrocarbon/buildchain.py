import sys
import os
import numpy as np

import build_HC_chain as bhc

#############################
## Number of repeating units
#############################

n=10

#################
## File to read
#################

outputfile='HC_chain_n'+str(n)+'.xyz'

#################
## Call Function
#################

bhc.build_full_HC_chain(outputfile,n)




