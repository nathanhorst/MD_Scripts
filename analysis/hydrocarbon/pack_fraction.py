# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 09:13:36 2016

@author: waltmann
"""

import numpy as np 
import V_distance as vd
import util as u



    
def polymer_volume(inputfile):
    distance_scale=1.0/3.96
    s=vd.type_pos_matrix(inputfile,'S')
    c2=vd.type_pos_matrix(inputfile,'CH2')
    c3=vd.type_pos_matrix(inputfile,'CH3')
    sigmas=np.array([4.45,3.96,3.76])
    volume=0
    volume+=(len(s)-1)* (sigmas[0]*distance_scale)**3
    volume+=(len(c2)-1)* (sigmas[1]*distance_scale)**3
    volume+=(len(c3)-1)* (sigmas[2]*distance_scale)**3
    return volume


####assumes nps are scaled where r=2 or some integer factor of that, also all
####np have the same number Au    
def np_volume(inputfile):
    au=vd.type_pos_matrix(inputfile,'Au')
    v=vd.type_pos_matrix(inputfile,'V')
    numau=(len(au)-1)/(len(v)-1)
    volume=0
    for x in range(1,len(v)):        
        for i in range(1+(numau*x-1),1+numau*x):
            r=u.part_distance(au[i],v[x],inputfile)
            if(r-float(int(r))==0.0):
                volume+=0.5*r**3
            i=-1
            print i
    return volume


###r is an array of all the different Rs probably just 2 times a scaling factor, 
####numAu is an array with corresponding nuber of Aus and 
### numnp is the corresponding number of NPs for that type           
def np_volume_known(r,numAu,numnp):
    volume=0.0
    for i in range(0,len(r)):
        volume+=numnp[i]*.5*r[i]**3
    return volume
        
def packing_fraction(inputfile):
    l=u.read_box_dimensions(inputfile)
    print l
    box=l[0]*l[1]*l[2]
    vol=(np.pi/6)*polymer_volume(inputfile)+np_volume_known([2.0],[122],[1])
    return vol/box
print packing_fraction('atoms.dump.0001500000.xml')
print packing_fraction('atoms.dump.0003480000.xml')
