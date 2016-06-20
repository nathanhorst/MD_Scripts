# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:09:51 2016

@author: waltmann
"""
import os
import sys

import numpy as np
from numpy import linalg as LA
##inputfile is xml



def vector_length(x,y,z):
    c=np.sqrt(x**2+y**2+z**2)
    return c
    
##assumes all chains are the same length, start with a sulfur, and the chains are
##all together at the end of the file in terms of particle number
    
 ###n can be any number 1-(length of chain(including sulfur)-1)   
def distance_to_nth_carbon(inputfile,nth):
    
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    active=0
    chnum=0
    cnum=0
    numAu=0
    box_size=0
    for line in data1:
        s=line.split()
        if(s[0]=='Au' or s[0]=='V'):
            numAu+=1
        elif(s[0]=='S'):
            active+=1
            chnum+=1
        elif(active==1):
            cnum+=1
        elif(s[0][:4]=='<box'):
            box_size=float(s[2][4:-1])
    count=0
    #print numAu
    #print cnum
    #print chnum
    pos=0
    #n=np.array([[0.0,0.0]])
    carnum=[]
    distance=[]
    spos=[0,0,0]   
    data2=f1data.splitlines()
    for line in data2:
        s=line.split()
        if(pos==1):
            count+=1
        if(s[0]=='<position'):
            pos+=1
        elif(s[0]=='</position>'):
            pos+=1
        elif(count>numAu and pos==1):
            if((count-numAu)%(cnum+1)!=1):
                carnum.append(int((count-numAu-1)%(cnum+1)))                
                x=abs(float(s[0])-spos[0])
                if(x>(box_size/2)-1):
                    x=box_size-x
                y=abs(float(s[1])-spos[1])
                if(y>(box_size/2)-1):
                    y=box_size-y
                z=abs(float(s[2])-spos[2])
                if(z>(box_size/2)-1):
                    z=box_size-z
                a=vector_length(x,y,z)
                distance.append(a)
            else:
                spos=[float(s[0]),float(s[1]), float(s[2])]
    distr=[]
    lc=0
    for u in range(0,len(carnum)):
        if(carnum[u]==nth):
            distr.append(distance[u])
            #print(str(distance[u]))
            lc+=1
            #print('entry')
        
    #print(carnum)
    #print(distance)
  #  n=np.array([[carnum[0]],distance[0]])
   # for i in range(1,len(carnum)):
    #    n=np.append(n,[carnum[i],distance[i]], axis=0)
    #print(distance[34])
    #return carnum,distance
    #print carnum
    #print distance
    distr.sort()
    #print distr
    #print len(carnum)    
    return distr        
#print(distance_to_nth_carbon('atoms.dump.0006475000.xml',12))    

 #n=np.append(n,[(count-numAu)%(cnum+1),vector_length(float(s[0])-spos[0],float(s[1])-spos[1], float(s[2])-spos[2])],axis=0)

def dist_average_of_files(nfiles,file1,file2, nth):
    avefile=0.0
    total=0.0
    file1=file1[11:-4]
    file2=file2[11:-4]
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        listed=distance_to_nth_carbon(toopen,nth)
        for i in range(0,len(listed)):
            total+=listed[i]
        avefile+=total/len(listed)
        total=0.0
    avefile=avefile/nfiles
    return avefile
#print(dist_average_of_files(16,'atoms.dump.0004600000.xml','atoms.dump.0004725000.xml', 12))