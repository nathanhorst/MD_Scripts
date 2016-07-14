# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:09:51 2016

@author: waltmann
"""
import os
import sys

import numpy as np
import V_distance as v
import util as u
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
    pos=0
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
            lc+=1
    distr.sort()  
    return distr          

def new_distance_to_nth_carbon(inputfile,nth):
    typs=['S','O','OH']
    x=nth_carbon_pos_matrix(inputfile,nth,typs)
    y=nth_carbon_pos_matrix(inputfile,0,typs)
    total=0
    for i in range(1,len(y)):
        xn=np.array([float(x[i][0]),float(x[i][1]),float(x[i][2])])
        yn=np.array([float(y[i][0]),float(y[i][1]),float(y[i][2])])
        total+=u.part_distance(xn,yn,inputfile)
    return float(total)/float(len(y)-1)
    
def new_distance_to_nth_carbon_list(inputfile,nth):
    typs=['S','O','OH']
    list=[]
    x=nth_carbon_pos_matrix(inputfile,nth,typs)
    y=nth_carbon_pos_matrix(inputfile,0,typs)
    for i in range(1,len(y)):
        xn=np.array([float(x[i][0]),float(x[i][1]),float(x[i][2])])
        yn=np.array([float(y[i][0]),float(y[i][1]),float(y[i][2])])
        list.append(u.part_distance(xn,yn,inputfile))
    return list
    
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
    
def dist_v_timestep(nfiles,file1,file2, nth):
    file1=file1[11:-4]
    file2=file2[11:-4]
    out=np.array([[0.0,0.0]])
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        x=new_distance_to_nth_carbon(toopen,nth)
        out=np.append(out,[[file0,x]],axis=0)
    return out

######typs is the types that constitute the chain s,ch2,ch3 head,repeat, end
def chain_pos_matrix(inputfile,typs):
    chain_pos=np.array([['type',0.0,0.0,0.0]])
    s=v.type_pos_matrix(inputfile,typs[0])
    #print s
    c2=v.type_pos_matrix(inputfile,typs[1])
    c3=v.type_pos_matrix(inputfile,typs[2])
    scount=0
    c2count=0
    c3count=0
    for n in range(1,len(s)):
        chain_pos=np.append(chain_pos,[[typs[0],s[n][0],s[n][1],s[n][2]]],axis=0)
        scount+=1
        l=(len(c2)-1)/(len(s)-1)
        for i in range(0,l):
            chain_pos=np.append(chain_pos,[[typs[1],c2[i+(n-1)*l][0],c2[i+(n-1)*l][1],c2[i+(n-1)*l][2]]],axis=0)
            c2count+=1
        chain_pos=np.append(chain_pos,[[typs[2],c3[n][0],c3[n][1],c3[n][2]]],axis=0)
        c3count+=1
    return chain_pos
    

def nth_carbon_pos_matrix(inputfile,n,typs):
    chain=chain_pos_matrix(inputfile,typs)
    n_mat=np.array([[0.0,0.0,0.0]])
    d=2
    while(chain[d][0]!=typs[0]):
        d+=1
    length=d-1
    d=n+1
    while(d<len(chain)):
        n_mat=np.append(n_mat,[[chain[d][1],chain[d][2],chain[d][3]]],axis=0)
        d+=length
    return n_mat

