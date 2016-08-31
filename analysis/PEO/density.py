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
import rotate_vector as rv
#import nanoparticle_core as npc
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
    
def closest_matrix_from_V(vpos_mat, spos_mat,c2pos_mat,c3pos_mat, inputfile):
    close_array=np.array([1000])
    v=vpos_mat[1:]
    c2=c2pos_mat[1:]
    c3=c3pos_mat[1:]
    s=spos_mat[1:]
    n=len(c2)/(len(s))
    f=len(s)/len(v)
    for i in range(0,(len(v))):
         a=np.array([float(v[i][0]),float(v[i][1]),float(v[i][2])])
         for su in range(0,f):
                 b=np.array([float(s[su+i*f][0]),float(s[su+i*f][1]),float(s[su+i*f][2])])
                 close_array=np.append(close_array,[u.length(np.subtract(a,b))],axis=0)
                 for x in range(0,n):
                     b=np.array([float(c2[x+n*su+n*f*i][0]),float(c2[x+n*su+n*f*i][1]),float(c2[x+n*su+n*f*i][2])])
                     close_array=np.append(close_array,[u.length(np.subtract(a,b))],axis=0)
                 b=np.array([float(c3[su+i*f][0]),float(c3[su+i*f][1]),float(c3[su+i*f][2])])
                 close_array=np.append(close_array,[u.length(np.subtract(a,b))],axis=0)
    return close_array

def position_matrix_from_V(vpos_mat,Aupos_mat, spos_mat,c2pos_mat,c3pos_mat, inputfile):
    close_array=np.array([[100,100,100]])
    v=vpos_mat[1:]
    au=Aupos_mat[1:]
    c2=c2pos_mat[1:]
    c3=c3pos_mat[1:]
    s=spos_mat[1:]
    thetax=0
    thetay=0
    thetaz=0
    correct=[0,0,2]
    n=len(c2)/(len(s))
    f=len(s)/len(v)
    for i in range(0,(len(v))):
         a=np.array([float(v[i][0]),float(v[i][1]),float(v[i][2])])
         for d in range(i*122,(i+1)*122):
                 e=np.array([float(au[d][0]),float(au[d][1]),float(au[d][2])])
                 if(abs(u.length(np.subtract(a,e))-2)<.00001):
                     thetaz=rv.angle_difference(correct[1],correct[0],e[1],e[0])
                     thetay=rv.angle_difference(correct[2],correct[0],e[2],e[0])
                     thetax=rv.angle_difference(correct[2],correct[1],e[2],e[1])
                     d=-1
         for su in range(0,f):
             b=np.array([float(s[su+i*f][0]),float(s[su+i*f][1]),float(s[su+i*f][2])])
             close_array=np.append(close_array,[np.subtract(rv.move_point(thetax,thetay,thetaz,b),a)],axis=0)
             for x in range(0,n):
                 b=np.array([float(c2[x+n*su+n*f*i][0]),float(c2[x+n*su+n*f*i][1]),float(c2[x+n*su+n*f*i][2])])
                 close_array=np.append(close_array,[np.subtract(rv.move_point(thetax,thetay,thetaz,b),a)],axis=0)
             b=np.array([float(c3[su+i*f][0]),float(c3[su+i*f][1]),float(c3[su+i*f][2])])
             close_array=np.append(close_array,[np.subtract(rv.move_point(thetax,thetay,thetaz,b),a)],axis=0)
    return close_array

def density_from_V(inputfile,r_step,r_limit,volume_divison):
    v=vd.type_pos_matrix(inputfile,'V')
    s=vd.type_pos_matrix(inputfile,'S')
    c2=vd.type_pos_matrix(inputfile,'O') 
    c3=vd.type_pos_matrix(inputfile,'OH')
    r_min=0
    r_max=r_step
    profile=np.array([[-1.0,0.0]])
    closest=closest_matrix_from_V(v,s,c2,c3,inputfile)
    while(r_max<r_limit):
        if(volume_divison):
            profile=np.append(profile,[[r_min,one_density(closest,r_min,r_max)]],axis=0)
        else:
            profile=np.append(profile,[[r_min,one_density_no_volume(closest,r_min,r_max)]],axis=0)
        r_min=r_max
        r_max+=r_step
    return profile

def density_cubeO_from_V(inputfile,r_step,r_limit,volume_divison):
    v=vd.type_pos_matrix(inputfile,'V')
    s=vd.type_pos_matrix(inputfile,'S')
    c2=vd.type_pos_matrix(inputfile,'O') 
    c3=vd.type_pos_matrix(inputfile,'OH')
    au=vd.type_pos_matrix(inputfile,'Au')
    r_min=0
    r_max=r_step
    profile=np.array([[-1.0,0.0]])
    closest=position_matrix_from_V(v,au,s,c2,c3, inputfile)
    #print closest
    #print 'matrix calculated'
    while(r_max<r_limit):
        if(volume_divison):
            profile=np.append(profile,[[r_min,one_cubeO_density(closest,r_min,r_max)]],axis=0)
        else:
            profile=np.append(profile,[[r_min,one_cubeO_density_no_volume(closest,r_min,r_max)]],axis=0)
        r_min=r_max
        r_max+=r_step
    total=0
    for i in range(0,len(profile)):
        total+=float(profile[i][1])
    #print total
    return profile
    
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

def one_cubeO_density_no_volume(close_array,r_min,r_max):
    particles=0.0
    mi=npc.NanoFcc(2.5,r_min/2.0)
    ma=npc.NanoFcc(2.5,r_max/2.0)
    for i in range(1,len(close_array)):
        if(ma.check_point(close_array[i]) and not mi.check_point(close_array[i])):
            particles+=1
    del mi
    del ma
    return particles
        
def one_cubeO_density(close_array,r_min,r_max):
    particles=0.0
    mi=npc.NanoFcc(2.5,r_min/2.0)
    ma=npc.NanoFcc(2.5,r_max/2.0)
    for i in range(1,len(close_array)):
        if(ma.check_point(close_array[i]) and not mi.check_point(close_array[i])):
            particles+=1
    volume=0.5* (r_max**3-r_min**3)
    x=float(particles)/volume
    del mi 
    del ma
    return x   
        
distance_scale=1.0
#u.write_array(density('atoms.dump.0000799968.xml', distance_scale/5, 30*distance_scale,True),'dense.txt')

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
