# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 08:52:46 2016

@author: waltmann
"""
from __future__ import division
import os
import sys

import numpy as np
from numpy import linalg as LA

###############################
##defines arctan on the entrire interval from (0,2pi) in order be used in angle
##difference function
###
#### -3pi is just an error code corresponding to (0,0) inputs
###################

def arc_tan_360(opposite,adjacent):
    if(adjacent==0):
        if(opposite==0):
            return -3*np.pi
        elif(opposite>0):
            return np.pi/2
        else:
            return 3*np.pi/2
    elif(adjacent<0):
        return np.arctan(opposite/adjacent)+np.pi
    elif(opposite<0):
        return np.arctan(opposite/adjacent)+2*np.pi
    else:
        return np.arctan(opposite/adjacent)
    
    
##############
#arctan(opp1/ad1)-arctan(opp2/ad2)
#########

def angle_difference(opp1,ad1,opp2,ad2):
    theta1=arc_tan_360(opp1,ad1)
    theta2=arc_tan_360(opp2,ad2)
    if(theta1==-3*np.pi):
        return 0
    elif(theta2==-3*np.pi):
        return 0
    else:
        return theta1-theta2
    
    

######################
##aligns the poltymer chain to v and then returns the polymer chain as a list
########################

    
"""
calculates the shift needed to line up the polymer with the certain vector for 
all 3 2d planes. Shifts in other planes will throw off a certain plane slightly
so these shifts are done in a loop until some tolerance is met
""" 
def align_vector(readpoly,v):
    fin1 = open(readpoly,'r')
    #fout = open('testchain.xyz','w')
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[2:]
    P=np.array([[0.0,0.0,0.0]])
    count=0          
    for line in data1:
        d=line.split()
        pos= [[float(d[1]),float(d[2]),float(d[3])]]
        P=np.append(P,pos,axis=0)
        count+=1;
        #fout.write(("%s %f %f %f\n")%(d[0],pos[0],pos[1],pos[2]))
    first=P[1]
    last=P[3]
    pv=[float(last[0])-float(first[0]),float(last[1])-float(first[1]),float(last[2])-float(first[2])]
    #print pv    
    pv=pv/(np.sqrt(pv[0]**2+pv[1]**2+pv[2]**2))
    v=v/(np.sqrt(v[0]**2+v[1]**2+v[2]**2))
    #print pv
    tolerance=.0001
    inside=0
    while(inside==0):    
        thetaz=angle_difference(v[1],v[0],pv[1],pv[0])
        thetaz=-thetaz
        Rz=[[np.cos(thetaz),-1*np.sin(thetaz),0],[np.sin(thetaz),np.cos(thetaz),0],[0,0,1]]
        P=np.dot(P,Rz)
        pv=np.dot(pv,Rz)
        #print('theta z is ' + str(thetaz*180/np.pi))
        #print(pv)
        thetay=angle_difference(v[2],v[0],pv[2],pv[0])
        Ry=[[np.cos(thetay),0,np.sin(thetay)],[0,1,0],[-1*np.sin(thetay),0,np.cos(thetay)]]
        P=np.dot(P,Ry)
        pv=np.dot(pv,Ry)
        #print('theta y is ' + str(thetay*180/np.pi))
        #print(pv)
        thetax=angle_difference(v[2],v[1],pv[2],pv[1])
        thetax=-thetax
        Rx=[[1,0,0],[0,np.cos(thetax),-1*np.sin(thetax)],[0,np.sin(thetax),np.cos(thetax)]]
        P=np.dot(P,Rx)
        pv=np.dot(pv,Rx)
        #print('theta x is ' + str(thetax*180/np.pi))
        #print(pv)
        if(np.abs(pv[0]-v[0])<tolerance and np.abs(pv[1]-v[1])<tolerance and np.abs(pv[2]-v[2])<tolerance):
            inside=1
    #print(P)
    #fout.write(str(count)+'\n\n')
    #for i in range(1,count+1):
        #fout.write(('C %s %s %s\n')%(P[i][0],P[i][1],P[i][2]))
    return P
    