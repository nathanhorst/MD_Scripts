# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 15:16:44 2016

@author: waltmann
"""

import os
import sys

import numpy as np
from numpy import linalg as LA

import dihedral_angle as da

import torsional as t

########
##n  is the dihedral number and can only range from one to length of chain - 3
######

#####
##l is the length of the chain including end groups
#########


def all_nth_dihedral(inputfile,n, l):
    dhs=t.all_dihedral_angles(inputfile)
    x=np.array([0])
    for v in range(0,len(dhs)):
        if((v%(l-3))==n-1):
            x=np.append(x,[dhs[v]],axis=0)
    return x
def all_nth_dihedral_lat(inputfile,n, l,totalnp,nptocount):
    dhs=t.some_dihedral_angles(inputfile,totalnp,nptocount)
    x=np.array([0])
    for v in range(0,len(dhs)):
        if((v%(l-3))==n-1):
            x=np.append(x,[dhs[v]],axis=0)
    return x    
    
def trans_nth_dihedral(inputfile,n, l):
    dhs=all_nth_dihedral(inputfile,n, l)
    trans=0.0
    for i in range(0,len(dhs)):
        if abs(dhs[i])<(np.pi/3):
            trans+=1
    return trans/len(dhs)
    
def trans_nth_dihedral_lat(inputfile,n, l,totalnp,nptocount):
    dhs=all_nth_dihedral_lat(inputfile,n, l,totalnp,nptocount)
    trans=0.0
    for i in range(0,len(dhs)):
        if abs(dhs[i])<(np.pi/3):
            trans+=1
    return trans/len(dhs)

def gplus_nth_dihedral(inputfile,n, l):
    dhs=all_nth_dihedral(inputfile,n, l)
    gplus=0.0
    for i in range(0,len(dhs)):
        if dhs[i]>(np.pi/3):
            gplus+=1
    return gplus/len(dhs)
    
def gminus_nth_dihedral(inputfile,n, l):
    dhs=all_nth_dihedral(inputfile,n, l)
    gminus=0.0
    for i in range(0,len(dhs)):
        if dhs[i]<(-1*np.pi/3):
            gminus+=1
    return gminus/len(dhs)

#print(all_nth_dihedral())
def average_trans(inputfile):
    numtrans=0.0
    numgplus=0.0
    numgminus=0.0    
    trans_cut=1.0472
    d=t.all_dihedral_angles(inputfile)
    print 'processing'
    for x in range(1,len(d)):
        if(abs(d[x])<trans_cut):
            numtrans+=1;
        elif(d[x]<-1.0*trans_cut):
            numgminus+=1
        else:
            numgplus+=1
        #print(float(len(d)-1))   
        
    return numtrans/float(len(d)-1)
    
    
def average_trans_lat(inputfile,totalnp,nptocount):
    numtrans=0.0
    numgplus=0.0
    numgminus=0.0    
    trans_cut=1.0472
    d=t.some_dihedral_angles(inputfile,totalnp,nptocount)
    #print 'processing'
    for x in range(1,len(d)):
        if(abs(d[x])<trans_cut):
            numtrans+=1;
        elif(d[x]<-1.0*trans_cut):
            numgminus+=1
        else:
            numgplus+=1
        #print(float(len(d)-1))   
        
    return numtrans/float(len(d)-1)
#print average_trans_lat('atoms.dump.0001100000.xml',32,1)
#print('   '+str(average_trans('11chain.xml')))
#############
##this figures out the average percent trans or gauce over a number of files throughout a
####simultaion for either the nth or all dihedrals
#######################
######
##assumes files are output at regular time intervals
#####
###############
##nth is 0 for all dihedrals or the number of the dihedral
##### l is the length of the chains
##inclusing end groups (S/CH3)
###############
###############

##############

def trans_average_of_files(nfiles,file1,file2, nth, l):
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
        if(nth==0):
            for n in range(0,l-3):
                total+=trans_nth_dihedral(toopen,n+1,l)
            avefile+=total/(l-3)
            total=0.0
        else:
            avefile+=trans_nth_dihedral(toopen,nth,l)
    avefile=avefile/nfiles
    return avefile


#print(trans_average_of_files(19,'atoms.dump.0006541657.xml','atoms.dump.0006749990.xml',0, 13))

def trans_v_timestep(nfiles,file1,file2, nth, l):
    file1=file1[11:-4]
    file2=file2[11:-4]
    out=np.array([[0.0,0.0]])
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        #print(toopen)
        if(nth==0):
            out=np.append(out,[[file0,average_trans(toopen)]],axis=0)
        else:
            out=np.append(out,[[file0,trans_nth_dihedral(toopen,nth,l)]],axis=0)
        #out=np.delete(out,0,0)
    return out
    
def trans_v_timestep_lat(nfiles,file1,file2, nth, l,totalnp,nptocount):
    file1=file1[11:-4]
    file2=file2[11:-4]
    out=np.array([[0.0,0.0]])
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        #print(toopen)
        if(nth==0):
            out=np.append(out,[[file0,average_trans_lat(toopen,totalnp,nptocount)]],axis=0)
        else:
            out=np.append(out,[[file0,trans_nth_dihedral_lat(toopen,nth,l,totalnp,nptocount)]],axis=0)
        #out=np.delete(out,0,0)
    return out
"""   
x=trans_v_timestep(97,'atoms.dump.0000500000.xml','atoms.dump.0000552083.xml',0, 13)
for i in range(1,len(x)):
    print(str(x[i][0])+ ' '+ str(x[i][1]))
"""


def angle_v_timestep(nfiles,file1,file2, nth, l,):
    file1=file1[11:-4]
    file2=file2[11:-4]
    out=np.array([[0.0,0.0]])
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        #print(toopen)
        if(nth==0):
            out=np.append(out,[[file0,trans_nth_dihedral(toopen,nth,l)]],axis=0)
        else:
            out=np.append(out,[[file0,all_nth_dihedral(toopen,nth,l)]],axis=0)
        #out=np.delete(out,0,0)
    return out
    
"""    
x=angle_v_timestep(145,'atoms.dump.0002499984.xml','atoms.dump.0002552067.xml',1, 4)
numtrans=0
for i in range(1,len(x)):
    print(str(x[i][0])+ ' '+ str(x[i][1]))
    if abs(x[i][1])<np.pi/3:
        numtrans+=1
print(numtrans/float(len(x)-1))
"""



 
    
