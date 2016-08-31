# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 10:54:15 2016

@author: waltmann
"""
import sys

from numpy import linalg as LA
import numpy as np
import util as u

"""
returns  a matrix containing the posiitons of all the V particle in the order
they are found in the file.
"""

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
            #print v_position    
            return v_position



"""
this function should return an array of nearest neighbor distances for 
 each particle with the position of the distnace corresponding to the
 original position of the particle in the input array
"""
def n_neighbor_vector_list(v_pos_mat,inputfile):
    distance=[[0.0,0.0,0.0]]
    for i in range(1,len(v_pos_mat)):
        min=10000
        minvec=[0.0,0.0,0.0]
        for x in range(1,len(v_pos_mat)):
            if u.part_distance(v_pos_mat[i],v_pos_mat[x],inputfile)<min and i!=x:
                minvec=np.subtract(v_pos_mat[x],v_pos_mat[i])
                min=u.part_distance(v_pos_mat[i],v_pos_mat[x],inputfile)
        distance=np.append(distance,[minvec],axis=0)
    return distance
 
"""
returns a list of all the particle sneaet neighbor ditstnace in the order the  particles
are found in the V-pos-matrix/file
"""
def n_neighbor_distance_list(v_pos_mat,inputfile):
     distance=[]
     for i in range(1,len(v_pos_mat)):
        min=10000
        for x in range(1,len(v_pos_mat)):
            if u.part_distance(v_pos_mat[i],v_pos_mat[x],inputfile)<min and i!=x:
                min=u.part_distance(v_pos_mat[i],v_pos_mat[x],inputfile)
        distance.append(min)
     return distance
    
#print(n_neighbor_distance_list(v_pos_matrix('Lat108_25_Fcc.xml'),'Lat108_25_Fcc.xml'))   
"""
meant to be used for pair force as it only does the distance between the first two particles in the file/matrix
"""

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
        out=np.append(out,[[file0,u.part_distance(x[1],x[2],toopen)]],axis=0)
    return out
    
#u.write_array(dist_v_timestep(34,'atoms.dump.0001500000.xml','atoms.dump.0001520000.xml'),'distance.txt')
"""
V_pos_matrix ,but for whatever type you specify. typ is a string that is the same as who it is in the file
"""
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
"""
basically returns an array of whatever position in the file of all atoms the atoms of that specific type
are. useful for 2.0 rigid body declaration.
"""            
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

#print v_pos_matrix('initialparticle_shift.xml')
