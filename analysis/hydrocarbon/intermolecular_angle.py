# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 12:02:44 2016

@author: waltmann
"""


import numpy as np

import util as u
import V_distance as vd
import distance_nth_carbon as dnc

def angle_between(vec1,vec2):
    x=np.dot(vec1,vec2)/(u.length(vec1)*u.length(vec2))
    if(x>1.0):
	x=1.0
    #print u.length(vec1)
    #print u.length(vec2)
    theta=np.arccos(x)
    theta=np.degrees(theta)
    #print theta
    return theta

"""
returns a list of intermolecular angles from the v particles to each S to the nth carbon
List should be analyzed using some sort of distribution function in util or radial_distribution
"""
def inter_angle(inputfile,n):
    typs=['S','CH2','CH3']
    angle_list=np.array([-500])
    vector_list=np.array([[0.0,0.0,0.0]])
    s=dnc.nth_carbon_pos_matrix(inputfile,0,typs)
    v=vd.v_pos_matrix(inputfile)
    if (n!=0):
        c=dnc.nth_carbon_pos_matrix(inputfile,n,typs)
        for i in range(1,len(s)):
            sn=np.array([float(s[i][0]),float(s[i][1]),float(s[i][2])])
            cn=np.array([float(c[i][0]),float(c[i][1]),float(c[i][2])])
            g=(len(c))/(len(v)-1)
            d=1+ int(i-1)/g
            vn=u.part_vector(cn,sn,inputfile)
            vector_list=np.append(vector_list,[vn],axis=0)
        #print vector_list
        for x in range(1,len(vector_list)):
            for y in range(x+1,len(vector_list)):
                angle_list=np.append(angle_list,[angle_between(vector_list[x],vector_list[y])],axis=0)
    else:
        for m in range(1,13):
            c=dnc.nth_carbon_pos_matrix(inputfile,m,typs)
            for i in range(1,len(s)):
                sn=np.array([float(s[i][0]),float(s[i][1]),float(s[i][2])])
                cn=np.array([float(c[i][0]),float(c[i][1]),float(c[i][2])])
                g=(len(c))/(len(v)-1)
                d=1+ int(i-1)/g
                vn=np.subtract(cn,sn)
                vector_list=np.append(vector_list,[vn],axis=0)
        for x in range(1,len(vector_list)):
            for y in range(x+1,len(vector_list)):
                angle_list=np.append(angle_list,[angle_between(vector_list[x],vector_list[y])],axis=0)
    
    angle_list=angle_list[1:]
    print angle_list
    return angle_list

"""
TODO
"""
"""
def inter_angle_all_files(file_list,n):
    angle_list=[-500]    
    for i in range(1,len(file_list)):
        angle_list=np.append(angle_list,inter_angle('atoms.dump.''file_list[i],n),axis=0)
    return angle_list[1:]
"""
