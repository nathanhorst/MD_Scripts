# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 09:43:29 2016

@author: nathan
"""
import numpy as np

def center_lat_xyz(save='centered.xyz',read='base.xyz'):
 
    #Open input file for reading in xyz format
 
    fin=open(read,'r')
    fdata=fin.read()
    datain=fdata.splitlines()
    num=datain[0]
    data=datain[2:]
    
    vlow=np.array([0.0,0.0,0.0])
    vhigh=np.array([0.0,0.0,0.0])
    
    for line in data:
        s=line.split()
        x=float(s[1])
        y=float(s[2])
        z=float(s[3])
        if x <= vlow[0]:                          
            vlow[0]=x
        if y <= vlow[1]:                          
            vlow[1]=y
        if z <= vlow[2]:                          
            vlow[2]=z
        if x > vhigh[0]:                          
            vhigh[0]=x
        if y > vlow[1]:                          
            vhigh[1]=y
        if z > vhigh[2]:                          
            vhigh[2]=z
    print('vLOW')
    print vlow
    print('vHIGH')
    print vhigh
    
    vL=np.add(vhigh,np.abs(vlow))    
    
    print('vL')
    print vL
    
    shiftx=-vL[0]/2.-vlow[0]
    shifty=-vL[1]/2.-vlow[0]
    shiftz=-vL[2]/2.-vlow[0]
    #Open output file for writing in xyz format
    
    fout=open(save,'w')
    fout.write(('%s \n \n')%(num))
 
    for line in data:
        s=line.split()
        x=float(s[1])+shiftx
        y=float(s[2])+shifty
        z=float(s[3])+shiftz
        fout.write(("%s %f %f %f\n")%(s[0],x,y,z))

