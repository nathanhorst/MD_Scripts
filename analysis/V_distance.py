# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 10:54:15 2016

@author: waltmann
"""
import sys

from numpy import linalg as LA
import numpy as np
import nth_dihedral as nth

def v_pos_matrix(inputfile):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[4:]     
    
    tipdone=False
    posdone=False
    pos=False
    tip=False
    tipcount=0
    poscount=0
    vnums=np.array([-5])
    position= np.array([[0.0,0.0,0.0]])
    v_position= np.array([[0.0,0.0,0.0]])
    for line in data1:
        #print line
        if(tip):
            tipcount+=1
        if(pos):
            poscount+=1
        s=line.split()
        if(s[0]=='</type>'):
            tip=False
            tipdone=True
        elif(tip==True and s[0]=='V'):
            vnums=np.append(vnums,[tipcount],axis=0)
        elif(s[0][:5]=='<type'):
            tip=True
        elif(s[0]=='</position>'):
            pos=False 
            posdone=True
        elif(pos):
            vec=[float(s[0]),float(s[1]),float(s[2])]
            position=np.append(position,[vec], axis=0)
        elif(s[0][:9]=='<position'):
            pos=True
        elif(posdone and tipdone):
            for x in range(1,len(vnums)):
                v_position=np.append(v_position,[position[vnums[x]]], axis=0)
            return v_position



def vector_length(x,y,z):
    c=np.sqrt(x**2+y**2+z**2)
    return c
    
def length(vector):
    return vector_length(vector[0],vector[1],vector[2])    

def part_distance(part1,part2,lx,ly,lz):
    
    x=abs(float(part1[0])-float(part2[0]))
    if(x>(lx/2)-1):
        x=lx-x
    y=abs(float(part1[1])-float(part2[1]))
    if(y>(ly/2)-1):
        y=ly-y
    z=abs(float(part1[2])-float(part2[2]))
    if(z>(lz/2)-1):
        z=lz-z
    return vector_length(x,y,z)


def part_distance(part1,part2,inputfile):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    #print data1[3]
    s=data1[3].split()
    #print s
    lx=0
    ly=0
    lz=0
    for i in range(len(s)):
        if s[i][:2]=='lx':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            lx=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='ly':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            ly=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='lz':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            lz=float(s[i][count+1:b+count-len(s[i])])
    #print lx
    #print ly
    #print lz
    #print part1
    #print part2[0]
    x=abs(float(part1[0])-float(part2[0]))
    if(x>(lx/2)-1):
        x=lx-x
    y=abs(float(part1[1])-float(part2[1]))
    #print y
    if(y>(ly/2)-1):
        y=ly-y
    z=abs(float(part1[2])-float(part2[2]))
    if(z>(lz/2)-1):
        z=lz-z
    #print x
    #print y 
    #print z
    return vector_length(x,y,z)

###this function should return an array of nearest neighbor distances for 
## each particle with the position of the distnace corresponding to the
### original position of the particle in the input array

def n_neighbor_vector_list(v_pos_mat):
    distance=[[0.0,0.0,0.0]]
    for i in range(1,len(v_pos_mat)):
        min=10000
        minvec=[0.0,0.0,0.0]
        for x in range(1,len(v_pos_mat)):
            if part_distance(v_pos_mat[i],v_pos_mat[x],50,50,50)<min and i!=x:
                minvec=v_pos_mat[x]
        distance=np.append(distance,[minvec],axis=0)
    return distance
    
def dist_v_timestep(nfiles,file1,file2):
    file1=file1[11:-4]
    file2=file2[11:-4]
    out=np.array([[0.0,0.0]])
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        x=v_pos_matrix(toopen)
        out=np.append(out,[[file0,part_distance(x[1],x[2],55,55,55)]],axis=0)
    return out
    
#nth.write_array(dist_v_timestep#(100,'atoms.dump.0000000000.xml','atoms.dump.0000100000.xml'),'distance.txt')

def type_pos_matrix(inputfile,typ):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[4:]     
    
    tipdone=False
    posdone=False
    pos=False
    tip=False
    tipcount=0
    poscount=0
    vnums=np.array([-5])
    position= np.array([[0.0,0.0,0.0]])
    v_position= np.array([[0.0,0.0,0.0]])
    for line in data1:
        #print line
        if(tip):
            tipcount+=1
        if(pos):
            poscount+=1
        s=line.split()
        if(s[0]=='</type>'):
            tip=False
            tipdone=True
        elif(tip==True and s[0]==typ):
            vnums=np.append(vnums,[tipcount],axis=0)
        elif(s[0][:5]=='<type'):
            tip=True
        elif(s[0]=='</position>'):
            pos=False 
            posdone=True
        elif(pos):
            vec=[float(s[0]),float(s[1]),float(s[2])]
            position=np.append(position,[vec], axis=0)
        elif(s[0][:9]=='<position'):
            pos=True
        elif(posdone and tipdone):
            for x in range(1,len(vnums)):
                v_position=np.append(v_position,[position[vnums[x]]], axis=0)
            return v_position
def type_tags(inputfile,typ):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[4:]     
    tipdone=False
    tip=False
    tipcount=-1
    vnums=np.array([-5])
    for line in data1:
        #print line
        if(tip):
            tipcount+=1
        s=line.split()
        if(s[0]=='</type>'):
            tip=False
            tipdone=True
        elif(tip==True and s[0]==typ):
            vnums=np.append(vnums,[tipcount],axis=0)
        elif(s[0][:5]=='<type'):
            tip=True
        elif(tipdone):
            return vnums
