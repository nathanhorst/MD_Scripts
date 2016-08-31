# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:09:51 2016

@author: waltmann
"""
import os
import sys

import numpy as np
import V_distance as v
import util as u
##inputfile is xml



def vector_length(x,y,z):
    c=np.sqrt(x**2+y**2+z**2)
    return c
    
    
"""    
assumes all chains are the same length, start with a sulfur, and the chains are
all together at the end of the file in terms of particle number
n can be any number 1-(length of chain(including sulfur)-1) 
"""  
def distance_to_nth_carbon(inputfile,nth):
    
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    active=0
    chnum=0
    cnum=0
    numAu=0
    box_size=0
    for line in data1:
        s=line.split()
        if(s[0]=='Au' or s[0]=='V'):
            numAu+=1
        elif(s[0]=='S'):
            active+=1
            chnum+=1
        elif(active==1):
            cnum+=1
        elif(s[0][:4]=='<box'):
            box_size=float(s[2][4:-1])
    count=0
    pos=0
    carnum=[]
    distance=[]
    spos=[0,0,0]   
    data2=f1data.splitlines()
    for line in data2:
        s=line.split()
        if(pos==1):
            count+=1
        if(s[0]=='<position'):
            pos+=1
        elif(s[0]=='</position>'):
            pos+=1
        elif(count>numAu and pos==1):
            if((count-numAu)%(cnum+1)!=1):
                carnum.append(int((count-numAu-1)%(cnum+1)))                
                x=abs(float(s[0])-spos[0])
                if(x>(box_size/2)-1):
                    x=box_size-x
                y=abs(float(s[1])-spos[1])
                if(y>(box_size/2)-1):
                    y=box_size-y
                z=abs(float(s[2])-spos[2])
                if(z>(box_size/2)-1):
                    z=box_size-z
                a=vector_length(x,y,z)
                distance.append(a)
            else:
                spos=[float(s[0]),float(s[1]), float(s[2])]
    distr=[]
    lc=0
    for u in range(0,len(carnum)):
        if(carnum[u]==nth):
            distr.append(distance[u])
            lc+=1
    distr.sort()  
    return distr          

"""
does the distnance_nth_carbon calculation based on the type_pos_matrix function found in 
V_distance.
Actually calls nth_carbon_matrix, but data processing can be traced back to V_distance
"""
def new_distance_to_nth_carbon(inputfile,nth):
    typs=['S','CH2','CH3']
    x=nth_carbon_pos_matrix(inputfile,nth,typs)
    y=nth_carbon_pos_matrix(inputfile,0,typs)
    total=0
    for i in range(1,len(y)):
        xn=np.array([float(x[i][0]),float(x[i][1]),float(x[i][2])])
        yn=np.array([float(y[i][0]),float(y[i][1]),float(y[i][2])])
        total+=u.part_distance(xn,yn,inputfile)
    return float(total)/float(len(y)-1)


"""
does the same thing as new_distance_nth_carbon, but returns a list of all the distances.
can be graphed as a distribution using util or radial_distribution
"""    
def new_distance_to_nth_carbon_list(inputfile,nth):
    typs=['S','CH2','CH3']
    list=[]
    x=nth_carbon_pos_matrix(inputfile,nth,typs)
    y=nth_carbon_pos_matrix(inputfile,0,typs)
    for i in range(1,len(x)):
        xn=np.array([float(x[i][0]),float(x[i][1]),float(x[i][2])])
        yn=np.array([float(y[i][0]),float(y[i][1]),float(y[i][2])])
        list.append(u.part_distance(xn,yn,inputfile))
    return list

"""
n is the number of the longest carbon using the convention that S is 0
Meant to be used once to determine a nanoparticles hydrodynamic radius 
for initial placement in a lattice
"""    
def max_radius(inputfile,n):
    typs=['S','CH2','CH3']
    x=nth_carbon_pos_matrix(inputfile,n,typs)
    y=v.v_pos_matrix(inputfile)
    max=0
    for i in range(0,len(x)):
        d=u.part_distance(y[1],x[i],inputfile)
        if (d>max):
            max=d
    return max
  
"""
finds the average nth distance over some files. file1 should be lower in timestep.
Also assumes that all files were dumped periodically so it can 
compute what other files exist.
Good for constant runs not for interpolations
"""  
  
def dist_average_of_files(nfiles,file1,file2, nth):
    avefile=0.0
    total=0.0
    file1=file1[11:-4]
    file2=file2[11:-4]
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        listed=distance_to_nth_carbon(toopen,nth)
        for i in range(0,len(listed)):
            total+=listed[i]
        avefile+=total/len(listed)
        total=0.0
    avefile=avefile/nfiles
    return avefile

"""
This is good for interpolations when there is only one particle.
For lattices this will take very long. About 1 minute/particle/file
So for 5 files with 32 particles, over 160 minutes
Also not good for calling in a loop for all 
"""
    
def dist_v_timestep(nfiles,file1,file2, nth):
    file1=file1[11:-4]
    file2=file2[11:-4]
    out=np.array([[0.0,0.0]])
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        x=new_distance_to_nth_carbon(toopen,nth)
        out=np.append(out,[[file0,x]],axis=0)
    return out

"""
This is the function to use for lattices. 
totalnp specifies how many nanaoparticles are in the file and 
nptocount tells how many to actually analyze in order to speed up the process
Again, analyzing one particle in one file takes a little over a minute (assuming 13 long 80 sites)
length is the length of the chain including end groups and the S.
Needs to be used in conjunction with nth_distance_lat and n_files in util to get 
the same full matrix result as dist_v_timestep.
Good for running nth_distance_lat in loops overall distances, because this function only 
needs to be called once and it takes all the time. 
"""
    
def dist_v_timestep_lat_matrix(nfiles,file1,file2,length,totalnp,nptocount):
    typs=['S','CH2','CH3']    
    file1=file1[11:-4]
    file2=file2[11:-4]
    out=np.zeros((nfiles,nptocount*(length-1)*80))
    for c in range(0,nfiles):
        file0=(int(file1)+(int(file2)-int(file1))*c)
        toopen='atoms.dump.'
        for v in range(0,10-len(str(file0))):
            toopen=toopen + '0'
        toopen=toopen + str(file0) + '.xml'
        print('Calculating Carbon Distance for timestep '+str(file0)+ ' '+ str(c+1)+'/' +str(nfiles))
        ch=chain_pos_lat_matrix(toopen,typs,totalnp,nptocount)
        #print ch
        spos=[0,0,0]
        cpos=[0,0,0]
        count=0
        for i in range(1,len(ch)):
            if(ch[i][0]==typs[0]):
                spos=[float(ch[i][1]),float(ch[i][2]),float(ch[i][3])]
            else:
                cpos=[float(ch[i][1]),float(ch[i][2]),float(ch[i][3])]
                out[c][count]=u.part_distance(spos,cpos,toopen)
                count+=1
    return out
    
    
"""
dvtlm is the dist_v_timestep_lat_matirix
This is good to run in a loop to process all the nth carbons in the dvtlm,
because it doesn't process the files again so its much faster
"""
def nth_distance_lat(file_of_files,n,dvtlm,length):
    towrite=np.array([['timestep',0.0]])
    for i in range(1,len(file_of_files)):
        total=0
        count=0
        for x in range(dvtlm.shape[1]):
            if((x+1)%(length-1)==n%(length-1)):
                total+=dvtlm[i-1][x]
                #print x
                count+=1
        towrite=np.append(towrite,[[file_of_files[i],float(total)/float(count)]],axis=0)
    return towrite
        
"""
typs is the types that constitute the chain s,ch2,ch3 head,repeat, end
this uses the V_distance function type_pos matrix to read the parts of the chains out from a file
and then puts them in the correct order in the chain. 
"""
def chain_pos_matrix(inputfile,typs):
    chain_pos=np.array([['type',0.0,0.0,0.0]])
    s=v.type_pos_matrix(inputfile,typs[0])
    c2=v.type_pos_matrix(inputfile,typs[1])
    c3=v.type_pos_matrix(inputfile,typs[2])
    scount=0
    c2count=0
    c3count=0
    for n in range(1,len(s)):
        chain_pos=np.append(chain_pos,[[typs[0],s[n][0],s[n][1],s[n][2]]],axis=0)
        scount+=1
        l=(len(c2)-1)/(len(s)-1)
        for i in range(0,l):
            chain_pos=np.append(chain_pos,[[typs[1],c2[i+(n-1)*l][0],c2[i+(n-1)*l][1],c2[i+(n-1)*l][2]]],axis=0)
            c2count+=1
        chain_pos=np.append(chain_pos,[[typs[2],c3[n][0],c3[n][1],c3[n][2]]],axis=0)
        c3count+=1
    return chain_pos
    
"""
typs is the types that constitute the chain ['S','CH2','CH3'] [head,repeat, end]
This is like chain_pos_matrix, but only processes a certain number of nps
based on the arguments totalnp and nptocount, both of which must be integers
called by other _mat functions to speed up data processing by ignoring some of the data
"""
def chain_pos_lat_matrix(inputfile,typs,totalnp,nptocount):
    chain_pos=np.array([['type',0.0,0.0,0.0]])
    s=v.type_pos_matrix(inputfile,typs[0])
    #print 's'
    #print s
    c2=v.type_pos_matrix(inputfile,typs[1])
    #print 'c2'
    #print c2
    c3=v.type_pos_matrix(inputfile,typs[2])
    scount=0
    c2count=0
    c3count=0
    x=float(nptocount)/float(totalnp)
    for n in range(1,int((len(s)-1)*x+1)):
        chain_pos=np.append(chain_pos,[[typs[0],s[n][0],s[n][1],s[n][2]]],axis=0)
        scount+=1
        l=(len(c2)-1)/(len(s)-1)
        for i in range(0,l):
            chain_pos=np.append(chain_pos,[[typs[1],c2[1+i+(n-1)*l][0],c2[1+i+(n-1)*l][1],c2[1+i+(n-1)*l][2]]],axis=0)
            c2count+=1
        chain_pos=np.append(chain_pos,[[typs[2],c3[n][0],c3[n][1],c3[n][2]]],axis=0)
        c3count+=1
    #print 'chain'
    #print chain_pos
    return chain_pos    

"""
reads the through the chain pos matrix and picks out the n carbon. This and all other functions
in this file use the convention that S is carbon 0 and the count goes up as you head down the chain
"""
def nth_carbon_pos_matrix(inputfile,n,typs):
    chain=chain_pos_matrix(inputfile,typs)
    n_mat=np.array([[0.0,0.0,0.0]])
    d=2
    while(chain[d][0]!=typs[0]):
        d+=1
    length=d-1
    d=n+1
    while(d<len(chain)):
        n_mat=np.append(n_mat,[[chain[d][1],chain[d][2],chain[d][3]]],axis=0)
        d+=length
    return n_mat


#print new_distance_to_nth_carbon_list('scaled_np.xml',12) 

"""
this is meant to be called on a lattice using the final carbon as n in order
to determine if the lattice spaicng used was correct or not.
"""

def min_distance_between_nth_carbon(inputfile,n):
    typs=['S','CH2','CH3']
    chain=nth_carbon_pos_matrix(inputfile,n,typs)
    #print 'matrix'
    min = 10000
    for i in range(1,len(chain)):
        #print i
        for x in range(i+1,len(chain)):
            d=u.part_distance(chain[i],chain[x],inputfile)
            if (d<min):
                min=d
    return min
#print min_distance_between_nth_carbon('2nano_12.xml',12)
#print max_radius('initialparticle_shift.xml',12)
