# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 15:49:33 2016

@author: waltmann
"""



import os
import sys

import numpy as np
from numpy import linalg as LA
import distance_nth_carbon as dnc 
import V_distance as vd
import nth_dihedral as nd
###########
##r_limit is the highest radius to be check
#############

##############
###r_step is the step by which the radii are checked
#############

def density(inputfile,r_step,r_limit):
    typs=['S','CH2','CH3']
    c=dnc.chain_pos_matrix(inputfile,typs)
    a=vd.type_pos_matrix(inputfile,'Au')
    r_min=0
    r_max=r_step
    profile=np.array([[-1.0,0.0]])
    close=closest_matrix(c,a,inputfile)
  #  print close
    while(r_max<r_limit):
        profile=np.append(profile,[[r_min,one_density(close,r_min,r_max)]],axis=0)
        r_min=r_max
        r_max+=r_step
    return profile
    
def closest_matrix(ppos_mat,Aupos_mat,inputfile):
    close_array=np.array([1000])
    for i in range(1,len(ppos_mat)):
        a=np.array([float(ppos_mat[i][1]),float(ppos_mat[i][2]),float(ppos_mat[i][3])])
        closest_distance=100000
        for x in range(1,len(Aupos_mat)):
            c=np.array([float(Aupos_mat[x][0]),float(Aupos_mat[x][1]),float(Aupos_mat[x][2])])
            distance=vd.part_distance(a,c,inputfile)
            if(distance<closest_distance):
                closest_distance=distance
        close_array=np.append(close_array,[closest_distance],axis=0)
    return close_array
    


def one_density(close_array,r_min,r_max):
    particles=0.0
    for i in range(1,len(close_array)):
        if (close_array[i]>=r_min and close_array[i]<r_max):
            particles+=1
    volume=4.0*np.pi/3.0* (r_max**3-r_min**3)
    x=float(particles)/volume
    #print x
    return x
    
#distance_scale=1.0/3.96
#nd.write_array(density('2np.xml', distance_scale/2, 6*distance_scale),'dense.txt')

def inter_angle(inputfile,n):
    typs=['S','CH2','CH3']
    angle_list=np.array([-500])
    s=dnc.nth_carbon_pos_matrix(inputfile,0,typs)
    c=dnc.nth_carbon_pos_matrix(inputfile,n,typs)
    #print len(c)
    #print len(s)
    v=vd.v_pos_matrix(inputfile)
    for i in range(1,len(s)):
        sn=np.array([float(s[i][0]),float(s[i][1]),float(s[i][2])])
        cn=np.array([float(c[i][0]),float(c[i][1]),float(c[i][2])])
        g=(len(c))/(len(v)-1)
        d=1+ int(i-1)/g
        vn=np.array([float(v[d][0]),float(v[d][1]),float(v[d][2])])
        normal=np.subtract(sn,vn)
        inter=np.subtract(cn,sn)
        theta=np.arccos(np.dot(normal,inter)/(vd.length(normal)*vd.length(inter)))
        theta=np.degrees(theta)
        angle_list=np.append(angle_list,[theta],axis=0)
    angle_list=angle_list[1:]
    return angle_list

#print(inter_angle('2np.xml',2))        

def normal_distribution(arr):
    std=np.std(arr)
    mean=np.mean(arr)
    dist=np.array([[0.0,0.0]])
    d=1/np.sqrt(2*std*std*np.pi)
    for i in range(0,60): 
        x=(mean-3*std)+i*(std/10.0)
        a=np.exp(-1*((x-mean)**2)/(2*std**2))*d
        dist=np.append(dist,[[x,a]],axis=0)
    return dist    
    

    
def histogram(arr,boxes):    
    arr=np.sort(arr,axis=0,kind='mergesort')
    print arr
    hist=np.array([[0.0,0.0]])
    step= (arr[len(arr)-1]-arr[0])/boxes
    for i in range(0,boxes):
        min=arr[0]+i*step
        max=arr[0]+(i+1)*step
        hit=0
        for x in range(0,len(arr)):
            if (arr[x]>=min and arr[x]<max):
                hit+=1
        hist=np.append(hist,[[(min+max)/2,hit]],axis=0)
    print hist
    return hist
nd.write_array(normal_distribution(inter_angle('atoms.dump.0005000000.xml',3)),'dense.txt')