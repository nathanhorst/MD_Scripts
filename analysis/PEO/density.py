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
import util as u
###########
##r_limit is the highest radius to be check
#############

##############
###r_step is the step by which the radii are checked
#############

def density(inputfile,r_step,r_limit,volume_divison):
    typs=['S','O','OH']
    c=dnc.chain_pos_matrix(inputfile,typs)
    a=vd.type_pos_matrix(inputfile,'Au')
    r_min=0
    r_max=r_step
    profile=np.array([[-1.0,0.0]])
    close=closest_matrix(c,a,inputfile)
  #  print close
    while(r_max<r_limit):
        if(volume_divison):
            profile=np.append(profile,[[r_min,one_density(close,r_min,r_max)]],axis=0)
        else:
            profile=np.append(profile,[[r_min,one_density_no_volume(close,r_min,r_max)]],axis=0)
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
            distance=u.part_distance(a,c,inputfile)
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

def one_density_no_volume(close_array,r_min,r_max):
    particles=0.0
    for i in range(1,len(close_array)):
        if (close_array[i]>=r_min and close_array[i]<r_max):
            particles+=1
    return particles
        
#distance_scale=1.0/3.96
#u.write_array(density('2np.xml', distance_scale/10, 14*distance_scale,False),'dense.txt')

def inter_angle(inputfile,n):
    typs=['S','O','OH']
    angle_list=np.array([-500])
    s=dnc.nth_carbon_pos_matrix(inputfile,0,typs)
    v=vd.v_pos_matrix(inputfile)
    if (n!=0):
        c=dnc.nth_carbon_pos_matrix(inputfile,n,typs)
        for i in range(1,len(s)):
            sn=np.array([float(s[i][0]),float(s[i][1]),float(s[i][2])])
            cn=np.array([float(c[i][0]),float(c[i][1]),float(c[i][2])])
            g=(len(c))/(len(v)-1)
            d=1+ int(i-1)/g
            vn=np.array([float(v[d][0]),float(v[d][1]),float(v[d][2])])
            angle_list=np.append(angle_list,[angle_between(vn,sn,cn)],axis=0)
    else:
        for m in range(1,13):
            c=dnc.nth_carbon_pos_matrix(inputfile,m,typs)
            for i in range(1,len(s)):
                sn=np.array([float(s[i][0]),float(s[i][1]),float(s[i][2])])
                cn=np.array([float(c[i][0]),float(c[i][1]),float(c[i][2])])
                g=(len(c))/(len(v)-1)
                d=1+ int(i-1)/g
                vn=np.array([float(v[d][0]),float(v[d][1]),float(v[d][2])])
                angle_list=np.append(angle_list,[angle_between(vn,sn,cn)],axis=0)
        angle_list=angle_list[1:]
        return angle_list

      
def angle_between(v,s,c):
    normal=np.subtract(s,v)
    inter=np.subtract(c,s)
    theta=np.arccos(np.dot(normal,inter)/(vd.length(normal)*vd.length(inter)))
    theta=np.degrees(theta)
    return theta


#nd.write_array(normal_distribution(inter_angle('atoms.dump.0005000000.xml',3)),'dense.txt')
