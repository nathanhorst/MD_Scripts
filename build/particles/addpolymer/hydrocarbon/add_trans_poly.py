# -*- coding: utf-8 -*-

import sys
import os
from numpy import *
import numpy as np
import rotate_polymer as rp
################
##t
##################

##############3
##R is the distance to center of the square face
## assumin the AU201 np
############


#####################
##readpoly should have no S
###############################3

def graft_square_faces(save,readpoly,readnp,R,sulfur_distance):
    fin1 = open(readpoly,'r')
    fin2 = open(readnp,'r')
    fout = open(save,'w')
    
    
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[2:]
    
    f2data=fin2.read()
    data2=f2data.splitlines()
    data2=data2[2:]
    #####################
    #AuSmultiplies a unit vestor to control the AuS initial distance
    ## AuS^2+.25^2= actual AuS distance ^2
    ###############3
    AuS=1
    s_s=sulfur_distance/np.sqrt(2)
    for line in data2:
        s=line.split()
        fout.write(('%s %f %f %f\n')%(s[0],float(s[1]),float(s[2]),float(s[3])))
    
##once the positions of the middle gold atom in the square face is found, four
##chians will be started at the s_s offset away    
    
    for line in data2:   
        s=line.split()
        vec=np.sqrt((float(s[1])**2)+(float(s[2])**2)+(float(s[3])**2))
        if vec == R:
            for i in range(0,4):
                    if(float(s[1])!=0):
                        r2=s_s*(int(i<2)*(1-(2*i)))+float(s[2])
                        r3=s_s*int(i>1)*(5-(2*i))+float(s[3])
                        r1= AuS*float(s[1])/float(vec)+float(s[1]) 
                    elif(float(s[2])!=0):
                        r1=s_s*(int(i<2)*(1-(2*i)))+float(s[1])
                        r3=s_s*int(i>1)*(5-(2*i))+float(s[3])
                        r2=AuS*float(s[2])/float(vec)+float(s[2]) 
                    else:
                        r1=s_s*(int(i<2)*(1-(2*i)))+float(s[1])
                        r2=s_s*int(i>1)*(5-(2*i))+float(s[2])
                        r3=AuS*float(s[3])/float(vec)+float(s[3]) 
                    
                    v=[r1,r2,r3]
                    P=rp.align_vector(readpoly,v)
                    for x in range(P.shape[0]):
                        P[x]=np.add(v,P[x])
                        if(x==P.shape[0]-1):
                            fout.write(('CH3 %f %f %f\n')%(P[x][0],P[x][1],P[x][2]))
                        elif(x==0):
                            fout.write(('S %f %f %f\n')%(P[x][0],P[x][1],P[x][2]))
                        else:
                            fout.write(('CH2 %f %f %f\n')%(P[x][0],P[x][1],P[x][2]))        
                                        
    
    
