# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 15:48:05 2016

@author: waltmann
"""



import numpy as np
import distance_nth_carbon as dnc
import util as u




def angle(v1,v2,v3,inputfile):
    vector1=u.part_vector(v2,v1,inputfile)
    vector2=u.part_vector(v2,v3,inputfile)
    angle=np.arccos(np.dot(vector1,vector2)/(u.length(vector1)*u.length(vector2)))
    return np.degrees(angle)
#print angle([1.121697, 1.121697, 1.121697],[1.205307, 1.485720, 1.220558],[1.486018, 1.486035, 1.486035],'Au140HC.xml') 
    
def bond_angle_v_timestep_lat_matrix(nfiles,file1,file2,length,totalnp,nptocount,numsites):
    typs=['S','CH2','CH3']    
    file1=file1[11:-4]
    file2=file2[11:-4]
    out=np.zeros((nfiles,nptocount*(length-1)*numsites))
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        print('Calculating Bond Angle for timestep '+str(file0)+ ' '+ str(c+1)+'/' +str(nfiles))
        ch=dnc.chain_pos_lat_matrix(toopen,typs,totalnp,nptocount)
        #print ch
        beforepos=[0,0,0]
        cpos=[0,0,0]
        nextpos=[0,0,0]
        count=0
        valid=False
        for i in range(1,len(ch)):
            if(ch[i][0]==typs[0]):
                nextpos=[float(ch[i][1]),float(ch[i][2]),float(ch[i][3])]
                valid=False
            else:
                beforepos=cpos
                cpos=nextpos
                nextpos=[float(ch[i][1]),float(ch[i][2]),float(ch[i][3])]
                if valid:
                    out[c][count]=angle(beforepos,cpos,nextpos,toopen)
                    print angle(beforepos,cpos,nextpos,toopen)
                    count+=1
                valid=True
    return out

