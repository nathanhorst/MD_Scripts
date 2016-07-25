# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 14:50:30 2016

@author: waltmann
"""

import numpy as np





#########################
##this program will create a matrix for all the particles that contains the 
#fields [[type,posx,posy,posz,mass,body]]
#######################
def status_matrix(inputfile):
    mass_ar=np.array([0])
    body_ar=np.array([0])
    typ_ar=np.array(['x'])
    pos_ar=np.array([[0,0,0]])
    typ=False
    typDone=False
    mass=False
    massDone=False
    pos=False
    posDone=False
    body=False
    bodyDone=False
    count=0
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[4:]  
    for line in data1:
       if(mass):
           count+=1
       if(typ):
            count+=1
       if(body):
           count+=1
       if(pos):
           count+=1
       s=line.split()
       if(s[0]=='</mass>'):
            mass=False
            massDone=True
            count=0
       elif(mass==True):
           mass_ar=np.append(mass_ar,[float(s[0])],axis=0)
       elif(s[0][:5]=='<mass'):
            mass=True
       elif(s[0]=='</type>'):
            typ=False 
            typDone=True
            count=0
       elif(typ):
            typ_ar=np.append(typ_ar,[s[0]], axis=0)
       elif(s[0][:5]=='<type'):
            typ=True
       elif(s[0]=='</body>'):
            body=False
            bodyDone=True
            count=0
       elif(body==True):
           body_ar=np.append(body_ar,[int(s[0])],axis=0)
       elif(s[0][:5]=='<body'):
            body=True
       elif(s[0]=='</position>'):
            pos=False
            posDone=True
            count=0
       elif(pos==True):
           pos_ar=np.append(pos_ar,[[float(s[0]),float(s[1]),float(s[2])]],axis=0)
       elif(s[0][:5]=='<posi'):
            pos=True
       elif(posDone and bodyDone and massDone and typDone):
           status_matrix=np.array([[typ_ar[1],float(pos_ar[1][0]),float(pos_ar[1][1]),float(pos_ar[1][2]),float(mass_ar[1]),int(body_ar[1])]])
           for i in range(2,len(body_ar)):
               y=np.array([[typ_ar[i],float(pos_ar[i][0]),float(pos_ar[i][1]),float(pos_ar[i][2]),float(mass_ar[i]),int(body_ar[i])]])
               status_matrix=np.append(status_matrix,y,axis=0)
           return status_matrix   

def status_matrix_shifted(s,shiftx,shifty,shiftz):
    ##[[type,posx,posy,posz,mass,body]]    
    v=np.array([['type',0.0,0.0,0.0,0.0,0.0]])    
    for i in range(len(s)):
        x=float(s[i][1])+shiftx
        y=float(s[i][2])+shifty
        z=float(s[i][3])+shiftz
        v=np.append(v,[[s[i][0],x,y,z,s[i][4],s[i][5]]],axis=0)
    return v[1:]

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


####[['type,posx,posy,posz']]
def all_pos_matrix_xyz(inputfile):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[2:]
    b=np.array([['type',0.0,0.0,0.0]])
    for line in data1:
        s=line.split()
        s[1]=float(s[1])
        s[2]=float(s[2])
        s[3]=float(s[3])
        b=np.append(b,[[s[0],s[1],s[2],s[3]]],axis=0)
    return b
    
            
def get_points(base_file):
    points=type_pos_matrix(base_file,'C')
    points=points[1:]
    return points

def get_points_xyz(base_file):
    points=all_pos_matrix_xyz(base_file)
    points=points[1:]
    return points    

def make_lattice(points,particle_file,newfile):
    a=np.array([['type',0.0,0.0,0.0]])
    x=status_matrix(particle_file)
    for p in points:
        v=status_matrix_shifted(x,float(p[1]),float(p[2]),float(p[3]))
        #print v
        #print 'fdd'
        for i in range(0,len(x)):
            a=np.append(a,[[v[i][0],v[i][1],v[i][2],v[i][3]]],axis=0)
    a=a[1:]
    fout=open(newfile,'w')
    fout.write(('%s \n \n')%(len(a)))
    
    for b in range(0,len(a)):
        x=float(a[b][1])
        y=float(a[b][2])
        z=float(a[b][3])
        fout.write(("%s %f %f %f\n")%(a[b][0],x,y,z))

###particle1 is the bigger 'A' particle, particle2 is the smaller 'B' particle

def make_binary_lattice(points,particle1,particle2,newfile):
    a=np.array([['type',0.0,0.0,0.0]])
    x=status_matrix(particle1)
    y=status_matrix(particle2)
    A=points[0][0]
    for p in points:
        if (p[0]==A):
            v=status_matrix_shifted(x,float(p[1]),float(p[2]),float(p[3]))
            #print v
            for i in range(0,len(v)):
                a=np.append(a,[[v[i][0],v[i][1],v[i][2],v[i][3]]],axis=0)
                #print 'wereked'
        else:
            v=status_matrix_shifted(y,float(p[1]),float(p[2]),float(p[3]))
            #print v
            #print 'dsfsdfasdfdasf'
            for i in range(0,len(v)):
                a=np.append(a,[[v[i][0],v[i][1],v[i][2],v[i][3]]],axis=0)
    a=a[1:]
    fout=open(newfile,'w')
    fout.write(('%s \n \n')%(len(a)))
    
    for b in range(0,len(a)):
        x=float(a[b][1])
        y=float(a[b][2])
        z=float(a[b][3])
        fout.write(("%s %f %f %f\n")%(a[b][0],x,y,z))
        
def particle_crystal(basefile,particlefile,newfile):
    make_lattice(get_points_xyz(basefile),particlefile,newfile)
    
basefile='base.xyz'
particle_file='Au201shellHC.xml'
new_file='new.xyz'

particle_crystal(basefile,particle_file,new_file)
    
