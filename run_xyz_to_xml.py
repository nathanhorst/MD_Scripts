# -*- coding: utf-8 -*-

import os
import sys

import numpy as np
import polymerxyztoxml as px


####3#####
##box length x,y,z
############
l=[40,40,40]

##################
#length of polymer including S and CH3
################
n=13

save='shift_check.xml'

read='shifted_1.xyz'


px.polynp_conv(l,save,read,n) 
