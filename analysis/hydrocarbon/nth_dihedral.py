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

"""
n  is the dihedral number and can only range from one to length of chain - 3
n  is the dihedral number and can only range from one to length of chain - 3
"""
def all_nth_dihedral(inputfile,n, l):
    dhs=t.all_dihedral_angles(inputfile)
    x=np.array([0])
    for v in range(0,len(dhs)):
        if((v%(l-3))==n-1):
            x=np.append(x,[dhs[v]],axis=0)
    return x
    
    
"""
the version for lattices to eliminate some of the redundant data
called by other _lat functions
"""
def all_nth_dihedral_lat(inputfile,n, l,totalnp,nptocount):
    dhs=t.some_dihedral_angles(inputfile,totalnp,nptocount)
    x=np.array([0])
    for v in range(0,len(dhs)):
        if((v%(l-3))==n-1):
            x=np.append(x,[dhs[v]],axis=0)
    return x    

"""
percent trans of the nth_dihedral
"""
    
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

"""
percent gauche plus of the nth dihedral
"""

def gplus_nth_dihedral(inputfile,n, l):
    dhs=all_nth_dihedral(inputfile,n, l)
    gplus=0.0
    for i in range(0,len(dhs)):
        if dhs[i]>(np.pi/3):
            gplus+=1
    return gplus/len(dhs)

"""
percent gauche minus of the nth_dihedral
"""
    
def gminus_nth_dihedral(inputfile,n, l):
    dhs=all_nth_dihedral(inputfile,n, l)
    gminus=0.0
    for i in range(0,len(dhs)):
        if dhs[i]<(-1*np.pi/3):
            gminus+=1
    return gminus/len(dhs)

"""
percent trans over an entire file
"""
def average_trans(inputfile):
    numtrans=0.0
    numgplus=0.0
    numgminus=0.0    
    trans_cut=1.0472
    d=t.all_dihedral_angles(inputfile)
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


"""
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
"""


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


"""
this figures out the percent trans 
for each tiemstep file and then outputs an array like
[[0.0,0.0],[timestep1,percent trans for timestep 1],[ect.]]
first must be cut offf and automatically is in the util wirte_array function
If n=0, it will do all the dihedrals no matter which one they are
"""

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

"""
same as trans_v_timestep, but cuts off the deadweight using totalnp and nptocount
"""    
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
this function is different because it returns a matrix that contains rows containing 
all dihedral andgels of for one file that you told it to read using the lattice argument totalnp and np tocount
The columns are all the files you told it to count. becuase of this the timestep info must be regained using n_files
and trnas_v_timestep_combine. This will get a matrix like the ones in the other timestep functions. 
Examples in run_polymer_analysis
"""    
def trans_v_timestep_lat_matrix(nfiles,file1,file2,l,totalnp,nptocount):
    file1=file1[11:-4]
    file2=file2[11:-4]
    out=np.zeros((nfiles,nptocount*80))
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        #print(toopen)
        print('Calculating Dihedral Angles for timestep '+ str(file0)+' '+ str(c)+'/'+str(nfiles))
        x=t.some_dihedral_angles(toopen,totalnp,nptocount)
        for i in range(1,80*nptocount+1):
            out[c][i-1]=x[i]
        print file0
    return out
"""
this is like the function above but without the matrix specifications.
Still needs to be recombined using the other functions and 
is most efficient when you wanna figure out all the nths in a loop
"""

def trans_v_timestep_matrix(nfiles,file1,file2,l):
    file1=file1[11:-4]
    file2=file2[11:-4]
    out=np.zeros((nfiles,80))
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        #print(toopen)
        x=t.all_dihedral_angles(toopen)
        for i in range(1,81):
            out[c][i-1]=x[i]
        print file0
    return out
 
"""
the matrix is those trans_v_timestep matrixes and file_of_files is what you
get by calling n_files. The result is a matrix like in the normal trans_v_timestep()
"""
   
def trans_v_timestep_combine(file_of_files,nth,matrix,length):
    towrite=np.array([['timestep',0.0]])
    for i in range(1,len(file_of_files)):
        towrite=np.append(towrite,[[file_of_files[i],nth_trans(matrix,nth,length)[i]]],axis=0)
    return towrite


"""
msut be called with the same file parameters that you called the trans_v_timeste_x function with
in order to recombine correctly
"""
def n_files(nfiles,file1,file2):
    file1=file1[11:-4]
    file2=file2[11:-4]
    files=np.array(['filename'])
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        files=np.append(files,[file0],axis=0)
    return files

#print n_files(2,'atoms.dump.0000300000.xml','atoms.dump.0000500000.xml')  

"""  
matrix is the trans_v_timestep_lat_matrix
called by the combine functuion.
This is what actually does the analysis in the matrix versions
and figures the percent trans of the nth dihdral.
"""

def nth_trans(matrix,n,length):
    out=np.array([0.0])
    if (n==0):
        for y in range(matrix.shape[0]):
            total=0.0
            for i in range(matrix.shape[1]):
                if (np.abs(matrix[y][i])<=np.pi/3):
                    total+=1.0
            out=np.append(out,[total/matrix.shape[1]],axis=0)
    else:
        for y in range(matrix.shape[0]):
            total=0.0
            hits=0.0
            for i in range(matrix.shape[1]):
                if((i+1)%(length-3)==(n%(length-3))):
                    hits+=1.0
                    if (np.abs(matrix[y][i])<=np.pi/3):
                        total+=1.0
            out=np.append(out,[total/hits],axis=0)
    return out


"""
should return a list of the dihedral angles for the nth_dihedral
n=0 means all dihedrals.
Not 100% sure this is working
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
    


 
    
