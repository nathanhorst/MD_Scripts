# -*- coding: utf-8 -*-

import os
import sys

import numpy as np
import polymerxyztoxml as px


###################
##box length x,y,z
###################
L=[40,40,40]

##################
#length of polymer including S and CH3
################
n=22

save='Au201PEO.xml'
read='Au201PEO.xyz'

px.polynp_conv(L,save,read,n) 
