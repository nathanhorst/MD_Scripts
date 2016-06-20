# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 12:53:20 2016

@author: waltmann
"""

import os
import sys

import numpy as np
from numpy import linalg as LA

import dihedral_angle as da




def s2v(x,y,z):
    v=np.array([float(x),float(y),float(z)])
    return v


##returns a list of all calculated dihedral angles

def all_dihedral_angles(inputfile):
    
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    
    pos=False
    tip=False
    position= np.array([[0.0,0.0,0.0]])
    tipe=np.array(['a'])
    for line in data1:
        s=line.split()
        if(s[0]=='</type>'):
            tip=False
        elif(tip==True):
            tipe=np.append(tipe,[s[0]],axis=0)
        elif(s[0][:5]=='<type'):
            tip=True
        elif(s[0]=='</position>'):
            pos=False    
        elif(pos):
            vec=[s[0],s[1],s[2]]
            position=np.append(position,[vec], axis=0)
        elif(s[0]=='<position'):
            pos=True
    dihedral=False
    #print(position)
    #print(tipe)
    dhs=np.array([0.0])
    for i in range(1,len(tipe)):
        if(tipe[i]=='S'):
            dihedral=True
        if(dihedral==True):
            one=s2v(position[i][0],position[i][1],position[i][2])
            two=s2v(position[i+1][0],position[i+1][1],position[i+1][2])
            three=s2v(position[i+2][0],position[i+2][1],position[i+2][2])
            four=s2v(position[i+3][0],position[i+3][1],position[i+3][2])
            dhs=np.append(dhs,[da.dihedral_angle(one,two,three,four)],axis=0)
            if(tipe[i+3]=='CH3'):
                dihedral=False            
    return dhs  
            
#print(all_dihedral_angles('atoms.dump.0019599999.xml'))                       
