# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 12:39:28 2016

@author: nathan
"""
import os
import sys

import numpy as np

def build_full_HC_chain(save,n):
    
    fid = open(save,'w')
    
    CH2=[0.0,0.0,0.0]
    #S=[-1.58237,0.899165,0.0] #### ACTUAL S VALUES
    S=[-1.58237,-0.883031,0.0] #### ADJUSTED S VALUEs
    Spos=np.add(CH2,S)
    CH2vec1=[1.24946,0.883031,0.0]
    CH2vec2=[1.24946,-0.883031,0.0]
 
    
    fid.write(str(2+n))
    fid.write(('\n\nS %s %s %s\n')%(Spos[0],Spos[1],Spos[2]))
    current= Spos
    for  x in range(0,n):
        if(x==0):
            current=CH2
        elif(x%2==0):
            current=np.add(current,CH2vec1)
        else:
            current=np.add(current,CH2vec2)
        fid.write(('CH2 %s %s %s\n')%(current[0],current[1],current[2]))
    if(n%2 == 0):
        newcurrent=np.add(current,CH2vec1)
        fid.write(('CH3 %s %s %s\n')%(newcurrent[0],newcurrent[1],newcurrent[2]))
    else:
        newcurrent=np.add(current,CH2vec2)
        fid.write(('CH3 %s %s %s\n')%(newcurrent[0],newcurrent[1],newcurrent[2]))
    
