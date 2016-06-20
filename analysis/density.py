# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 15:49:33 2016

@author: waltmann
"""



import os
import sys

import numpy as np
from numpy import linalg as LA


###########
##r_limit is the highest radius to be check
#############

##############
###r_step is the step by which the radii are checked
#############

def density(inputfile, r_limit, r_step):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    
    vpos=np.array([0,0,0])
    ppos=np.array([[0,0,0]])    
    tip=False
    pos=False
    toV=0
    poscount=0;
    for line in data1:
        s=line.split()
        if(s[0]=='</type>'):
            tip=False
        elif(tip==True):
            toV+=1
            if(s[0]=='V'):
                tip=False
        elif(s[0][:5]=='<type'):
            tip=True
        elif(s[0]=='</position>'):
            pos=False    
        elif(pos):
            poscount+=1
            vec=np.array([s[0],s[1],s[2]])
            if(poscount==toV):
                vpos=vec
            else:
                ppos=np.append(ppos,[vec], axis=0)
        elif(s[0]=='<position'):
            pos=False
        final= np.array([[0,0]])
        r=r_step
    while(r<=r_limit):
        de=one_density(ppos,vpos,r)
        final=np.append(final,[[r,de]],axis=0)
        r+=r_step
    return final


def one_density(ppos,vpos,r):
    particles=1.0
    for i in range(1,len(ppos)):
        ppos=np.subtract(ppos[i],vpos)
        len_ppos=np.sqrt(ppos[i][0]**2+ppos[i][1]**2+ppos[i][2]**2)
        if(len_ppos>=r):
            particles+=1
    dense=particles/(4*np.pi/3*r**3)
    return dense
    
print(density('atoms.dump.0009599998.xml', 50, .5))