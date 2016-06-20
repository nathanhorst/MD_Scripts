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
    
    
def trans_nth_dihedral(inputfile,n, l):
    dhs=all_nth_dihedral(inputfile,n, l)
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

for i in range(1,10):
    print(str(trans_nth_dihedral('atoms2.dump.0004599999.xml',i,13)))



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


#print(trans_average_of_files(16,'atoms.dump.0004600000.xml','atoms.dump.0004725000.xml',0, 13))