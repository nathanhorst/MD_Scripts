# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:10:00 2016

@author: waltmann
"""

import V_distance as vd
import util as u
import numpy as np
import density as d

def distance_array_V(inputfile):
    x=vd.v_pos_matrix(inputfile)
    d=[0.0,0.0,0.0]
    v=0
    for i in range(1,len(x)):
        min=100
        temp=u.part_distance(d,x[i],inputfile)
        if(temp<min):
            v=i
            min=temp
    distance_array=np.array([0.0])
    for b in range(1,len(x)):
        if b!=i:
            c=u.part_distance(x[b],x[v],inputfile)
            distance_array=np.append(distance_array,[c],axis=0)
    distance_array=np.sort(distance_array)
    return distance_array

def distance_array_chains(inputfile):
    typs=['S','O','OH']
    s=vd.type_pos_matrix(inputfile,typs[0])
    c2=vd.type_pos_matrix(inputfile,typs[1])
    c3=vd.type_pos_matrix(inputfile,typs[2])
    v=vd.v_pos_matrix(inputfile)
    x=d.closest_matrix_from_V(v,s,c2,c3,inputfile)
    return x
        
def radial_distribution(matrix):
    v=matrix
    print v
    #print v
    upper_limit=v[len(v)-1]+1
    lower_limit=0
    step=int(upper_limit*10)
    #print v[1]
    #print lower_limit
    #print upper_limit
    #print step
    total=u.dirac_delta(v[1],lower_limit,upper_limit,step)
    for y in range(0,len(total)):
        total[y][1]=total[y][1]/((4/3)*np.pi*v[1]**3)
    for x in range(2,len(v)):
        c=u.dirac_delta(v[x],lower_limit,upper_limit,step)
        for i in range(0,len(c)):
            c[i][1]=c[i][1]/((4/3)*np.pi*v[x]**3)
            total[i][1]+=c[i][1]
    for p in range(0,len(total)):
        total[p][1]+=total[p][1]/float(len(v)-1)
        total[p][1]+=total[p][1]/float(len(v)-1)
    #print total
    return total

#u.write_array(radial_distribution(distance_array_chains('atoms.dump.0001500000.xml')),'dd.txt')        
