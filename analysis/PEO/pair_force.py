# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 16:37:49 2016

@author: waltmann
"""
import V_distance as vd
import numpy as np

def vector_length(x,y,z):
    c=np.sqrt(x**2+y**2+z**2)
    return c
    
def length(vector):
    return vector_length(vector[0],vector[1],vector[2]) 
    
    
    
 ####takes the xml file and creates a 2d array of masses and acccelerations in order   
    
def mass_acc_list(inputfile):
   fin1= open(inputfile,'r')
   f1data=fin1.read()
   data1=f1data.splitlines()
   data1=data1[4:]     
   massdone=False
   acceldone=False
   mass=False
   accel=False
   masscount=0
   accelcount=0
   mass_ar=np.array([0.0])
   accel_ar=np.array([0.0])
   ma= np.array([[0.0,0.0]])
   for line in data1:
       if(mass):
           masscount+=1
       if(accel):
            accelcount+=1
       s=line.split()
       if(s[0]=='</mass>'):
            mass=False
            massdone=True
       elif(mass==True):
           mass_ar=np.append(mass_ar,[float(s[0])],axis=0)
       elif(s[0][:5]=='<mass'):
            mass=True
       elif(s[0]=='</acceleration>'):
            accel=False 
            acceldone=True
       elif(accel):
            vec=[float(s[0]),float(s[1]),float(s[2])]
            accel_ar=np.append(accel_ar,[length(vec)], axis=0)
       elif(s[0][:13]=='<acceleration'):
            accel=True
       elif(massdone and acceldone):
            for x in range(1,masscount):
                ma=np.append(ma,[[mass_ar[x],accel_ar[x]]], axis=0)
            return ma
#print(np_seperate(2,mass_acc_list('atoms.dump.0007000000.xml')))
     
#uses a mass acceleration list to calculate the F1 vector for each np
         
def np_seperate(numbnp,arr):
    forces= np.array([[0,0.0]])
    t=((len(arr)-1)/numbnp)+1
    s=1
    for i in range(1,numbnp+1):
        fn=0
        for x in range(s,t):
            fn+=arr[x][0]*arr[x][1]
        forces=np.append(forces,[[i,fn]],axis=0)
        x=(len(arr)-1)/numbnp  
        s+=x
        t+=x
    return forces

def force_average_of_files(nfiles,file1,file2):
    total=0.0
    avefile=0.0
    file1=file1[11:-4]
    file2=file2[11:-4]
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        #print(toopen)
        x=np_seperate(2,mass_acc_list(toopen))
        y=np.subtract(x[2],x[1])
        z=vd.v_pos_matrix(toopen)
        w=np.subtract(z[2],z[1])
        ru=w/vd.length(w)
        f=np.array([y[1]*y[0],y[2]*y[0],y[3]*y[0]])
        total+=np.dot(f,ru)
    avefile=float(total)/float(nfiles)
    return avefile    
