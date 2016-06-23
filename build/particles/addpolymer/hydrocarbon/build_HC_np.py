# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 14:15:00 2016

@author: nathan
"""
import os
import sys

import numpy as np
from numpy import linalg as LA
import rotate_polymer as rp

def build_HC_chain_noS(save,n):
    
    fid = open(save,'w')
    
    CH2=[0.0,0.0,0.0]
    S=[-1.58237,0.899165,0.0]
    Spos=np.add(CH2,S)
    CH2vec1=[1.24946,0.883031,0.0]
    CH2vec2=[1.24946,-0.883031,0.0]
    
    fid.write(str(1+n)+'\n\n')
    current= Spos
    for  x in range(0,n):
        if(x==0):
            current=CH2
        if(x%2==0):
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

def graft_HC_chain(save,readpoly,readnp,R1,R2):
    
    fin1 = open(readpoly,'r')
    fin2 = open(readnp,'r')
    fout = open(save,'w')
    
    
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[2:]
    
    f2data=fin2.read()
    data2=f2data.splitlines()
    data2=data2[2:]
    
    sulfur=[-1.58237,0.899165,0.0]
    
    for line in data2:
        s=line.split()
        fout.write(('%s %f %f %f\n')%(s[0],float(s[1]),float(s[2]),float(s[3])))
        
    for line in data2:   
        s=line.split()
        vec=np.sqrt((float(s[1])**2)+(float(s[2])**2)+(float(s[3])**2))
        if vec == R1:
            r1=float(s[1])/float(vec)+float(s[1]) 
            r2=float(s[2])/float(vec)+float(s[2])
            r3=float(s[3])/float(vec)+float(s[3])
            #fout.write(("S %f %f %f\n")%(r1,r2,r3))
            v=[r1,r2,r3]
            P=rp.align_vector(readpoly,v)
            for i in range(P.shape[0]):
                P[i]=np.add(v,P[i])
                #if(i==P.shape[1]-1):
                if(i==P.shape[0]-1):
                    fout.write(('CH3 %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
                elif(i==0):
                    fout.write(('S %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
                else:
                    fout.write(('CH2 %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
        if vec == R2:
            r1=float(s[1])/float(vec)+float(s[1]) 
            r2=float(s[2])/float(vec)+float(s[2])
            r3=float(s[3])/float(vec)+float(s[3])
            #fout.write(("S %f %f %f\n")%(r1,r2,r3))
            v=[r1,r2,r3]
            P=rp.align_vector(readpoly,v)
            for i in range(P.shape[0]):
                P[i]=np.add(v,P[i])
                #if(i==P.shape[1]-1):
                if(i==P.shape[0]-1):
                    fout.write(('CH3 %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
                elif(i==0):
                    fout.write(('S %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))  
                else:
                    fout.write(('CH2 %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
                
            
            
        