# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 08:58:19 2016

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


#input the status matrix you created in status_matrix
#shifts the np further to the positive x direction, shifts should be negative to bring closer together
#going to return the modified status matrix     
#assumes two particles for now      
def particle_shift(status_matrix,shiftx,shifty,shiftz):
    numbnp=2    
    n=len(status_matrix)/numbnp
    which =1
    if status_matrix[0][2]>status_matrix[n][2]:
        which = 0
    for i in range(0,n):
        status_matrix[i+which*n][1]=float(status_matrix[i+which*n][1])+shiftx
        status_matrix[i+which*n][2]=float(status_matrix[i+which*n][2])+shifty
        status_matrix[i+which*n][3]=float(status_matrix[i+which*n][3])+shiftz
    return status_matrix

###############
#this program should take an atoms file and reinitialize it as a motionless xml 
#with the particles closer together
################
def atoms_shifted_xml(save,read,shiftx,shifty,shiftz):
    x=status_matrix(read)
    x=particle_shift(x,shiftx,shifty,shiftz)
    fid = open(save,'w')
    inputfile = open(read,'r')
    fdata=inputfile.read()
    
    data=fdata.splitlines()
    data_short=data[:4]
    data=data[4:]
    cou=1
    for line in data_short:
        if cou==3:
            fid.write('<configuration time_step="0">\n')
        else:   
            fid.write(line)
            fid.write('\n')
        cou+=1
    fid.write('<type>\n')
    for a in range(0,len(x)):
        fid.write(str(x[a][0])+'\n')
    fid.write('</type>\n')
    fid.write('<position>\n')
    for b in range(0,len(x)):
        fid.write(str(x[b][1])+' '+ str(x[b][2])+' '+str(x[b][3])+'\n')
    fid.write('</position>\n')
    fid.write('<mass>\n')
    for c in range(0,len(x)):
        fid.write(str(x[c][4])+'\n')
    fid.write('</mass>\n')
    fid.write('<body>\n')
    for c in range(0,len(x)):
        fid.write(str(x[c][5])+'\n')
    fid.write('</body>\n')      
    dihedral=False
    dihedralDone=False
    angle=False
    angleDone=False
    bond=False
    bondDone=False    
    for line in data:
        s=line.split()
        if (s[0][:9]=='<dihedral'):
            dihedral=True
            fid.write('<dihedral>'+ '\n')
        elif(dihedral):
            for g in range(0,len(s)):
                fid.write(s[g]+ " ")
            fid.write('\n')
            if(s[0][:10]=='</dihedral'):
                dihedral=False
                dihedralDone=True
        elif (s[0][:5]=='<bond'):
            bond=True
            fid.write('<bond>'+ '\n')
        elif(bond==True):
            for g in range(0,len(s)):
                fid.write(s[g]+ " ")
            fid.write('\n')
            if(s[0][:6]=='</bond'):
                bond=False
                bondDone=True
        elif (s[0][:6]=='<angle'):
            angle=True
            fid.write('<angle>'+ '\n')
        elif(angle):
            for g in range(0,len(s)):
                fid.write(s[g]+ " ")
            fid.write('\n')
            if(s[0][:7]=='</angle'):
                angle=False
                angleDone=True
        elif(angleDone and bondDone and dihedralDone):
            data_end=data[-2:]
            for line in data_end:
                fid.write(line+ '\n')
            return
    
    
#atoms_shifted_xml('try.xml','atoms.dump.0004600000.xml',-1,0,0)    
        
