# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 11:16:41 2016

@author: waltmann
"""

import os
import sys

import numpy as np
from numpy import linalg as LA

#### the p inputs are positions in 3 dimensional space [x,y,z]
##should go in order of closest to nanoparticle first

def dihedral_angle(p1,p2,p3,p4):
    b1 = -1.0*(p2 - p1)
    b2 = p3 - p2
    b3 = p4 - p3
    
    b1xb2= np.cross(b1,b2)
    b2xb3=np.cross(b2,b3)
    
    
    magb2= np.sqrt(b2[0]**2+b2[1]**2+b2[2]**2)
    nb2=b2/magb2
    
    dc=np.cross(b1xb2,b2xb3)
    cd=np.dot(b1xb2,b2xb3)
    
    ccd=np.dot(dc,nb2)
    
    x=np.arctan2(ccd,cd)
    
    return x
    
    
    
    