# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 10:09:06 2016

@author: waltmann
"""

from __future__ import division
kb=8.617e-5

import sys
import scipy
import numpy as np
sys.path.insert(0,'../')
import util as u 
import scipy.optimize as opt
import numpy as np



def Esplit(Hfile,numrs,boxesperR):
    with open("mylog.log") as f:
        data = f.readlines()
        data=data[2:]
    #print data[1].split()
    with open(Hfile) as dist:
        data2 = dist.readlines()
        data2=data2[2:]
    H=np.array([-3.0])
    search=np.array([-3.0])
    for row in data2:
        search=np.append(search,[float(row.split()[0])*3.96],axis=0)
    for x in range(0,len(data)):
        H=np.append(H,[float(data[x].split()[2])*0.002872],axis=0)
    array=np.array([0.0])
    index=0
    rnum=0
    H=H[1:]
    print len(H)
    search=search[1:]
    #print 'H array made'
    splitt=np.zeros((numrs*2,boxesperR+2))
    perwindow=len(H)/numrs
    cutoff=perwindow/2
    while(index<len(H)):
        if(float(search[index])!=-3.0*3.96 and index%perwindow>cutoff):
            array=np.append(array,[float(H[index])],axis=0)
            print 'enter'
        elif(float(search[index])==-3.0*3.96):
            print "Gonna make a histogram"
            array=array[1:]
            v=u.histogram(array,int(boxesperR))
            v=v[1:]
            Ni=0
            for c in range(0,len(v)):
                Ni+=v[c][1]
            print Ni
            for g in range(0,len(v)):
                splitt[2*rnum][g]=v[g][0]
                splitt[2*rnum+1][g]=v[g][1]
            rnum+=1
            print(len(array))
            array=np.array([0.0])
        index+=1
         
    return splitt

def LJEsplit(Hfile,numrs,boxesperR):
    with open("mylog.log") as f:
        data = f.readlines()
        data=data[2:]
    #print data[1].split()
    with open(Hfile) as dist:
        data2 = dist.readlines()
        data2=data2[2:]
    H=np.array([-3.0])
    search=np.array([-3.0])
    for row in data2:
        search=np.append(search,[float(row.split()[0])*3.96],axis=0)
    for x in range(0,len(data)):
        H=np.append(H,[float(data[x].split()[6])*0.002872],axis=0)
    array=np.array([0.0])
    index=0
    rnum=0
    H=H[1:]
    print len(H)
    search=search[1:]
    #print 'H array made'
    splitt=np.zeros((numrs*2,boxesperR+2))
    perwindow=len(H)/numrs
    cutoff=perwindow/2
    while(index<len(H)):
        if(float(search[index])!=-3.0*3.96 and index%perwindow>cutoff):
            array=np.append(array,[float(H[index])],axis=0)
            print 'enter'
        elif(float(search[index])==-3.0*3.96):
            print "Gonna make a histogram"
            array=array[1:]
            v=u.histogram(array,int(boxesperR))
            v=v[1:]
            Ni=0
            for c in range(0,len(v)):
                Ni+=v[c][1]
            print Ni
            for g in range(0,len(v)):
                splitt[2*rnum][g]=v[g][0]
                splitt[2*rnum+1][g]=v[g][1]
            rnum+=1
            print(len(array))
            array=np.array([0.0])
        index+=1
         
    return splitt
    
def BondEsplit(Hfile,numrs,boxesperR):
    with open("mylog.log") as f:
        data = f.readlines()
        data=data[2:]
    #print data[1].split()
    with open(Hfile) as dist:
        data2 = dist.readlines()
        data2=data2[2:]
    H=np.array([-3.0])
    search=np.array([-3.0])
    for row in data2:
        search=np.append(search,[float(row.split()[0])*3.96],axis=0)
    for x in range(0,len(data)):
        H=np.append(H,[float(data[x].split()[7])*0.002872],axis=0)
    array=np.array([0.0])
    index=0
    rnum=0
    H=H[1:]
    print len(H)
    search=search[1:]
    #print 'H array made'
    splitt=np.zeros((numrs*2,boxesperR+2))
    perwindow=len(H)/numrs
    cutoff=perwindow/2
    while(index<len(H)):
        if(float(search[index])!=-3.0*3.96 and index%perwindow>cutoff):
            array=np.append(array,[float(H[index])],axis=0)
            print 'enter'
        elif(float(search[index])==-3.0*3.96):
            print "Gonna make a histogram"
            array=array[1:]
            v=u.histogram(array,int(boxesperR))
            v=v[1:]
            Ni=0
            for c in range(0,len(v)):
                Ni+=v[c][1]
            print Ni
            for g in range(0,len(v)):
                splitt[2*rnum][g]=v[g][0]
                splitt[2*rnum+1][g]=v[g][1]
            rnum+=1
            print(len(array))
            array=np.array([0.0])
        index+=1
         
    return splitt
    
def E_not_split(Hfile,numrs,boxesperR):
    with open(Hfile) as f:
        data = f.readlines()
    with open(Hfile) as dist:
        data2 = dist.readlines()
    H=np.array([-3.0])
    
    for x in range(0,data.shape[0]):
        if(float(data2[x].split()[0])!=-3.0):
            H=np.append(H,[float(data[x].split()[2])*.002872],axis=0)
    H=H[1:]
    print H[1]
    v=u.histogram(H,int(boxesperR*numrs))
    v=v[1:]
    return v


def concat_h(file1,file2,file3):
      with open(file1) as f1:
        data1 = f1.readlines()
      H=np.array([-3.0])
      for row in data1:
        H=np.append(H,[float(row.split()[0])],axis=0)
      with open(file2) as f2:
          data2 = f2.readlines()
      for row in data2:
          H=np.append(H,[float(row.split()[0])],axis=0)
      H=H[1:]
      with open(file3,'w') as f3:
           for i in range(0,len(H)):
               f3.write("%f" % (H[i]))
               f3.write('\n')        

#concat_h("H.txt","fake.txt","H.txt")
def cut_lines_H(Hfile,lines,ofile):
    with open(Hfile) as f:
        data=f.readlines()
    H=np.array([-3.0])
    for row in data:
        H=np.append(H,[row.split()[0]],axis=0)
    H=H[lines+1:]
    with open(ofile,'w') as w:
        for i in range(0,len(H)):
            w.write(H[i])
            w.write('\n')


def cut_lines_H_end(Hfile,lines,ofile):
    with open(Hfile) as f:
        data=f.readlines()
    H=np.array([-3.0])
    for row in data:
        H=np.append(H,[row.split()[0]],axis=0)
    H=H[:len(H)-lines]
    with open(ofile,'w') as w:
        for i in range(0,len(H)):
            w.write(H[i])
            w.write('\n')
#cut_lines_H_end("18.txt",120012,"6.txt")         
    
def Hsplit(Hfile,numrs,boxesperR):
    with open(Hfile) as f:
        data = f.readlines()
    H=np.array([-3.0])
    
    for row in data:
        H=np.append(H,[float(row.split()[0])*3.96],axis=0)
    array=np.array([0.0])
    index=0
    rnum=0
    H=H[1:]
    #print 'H array made'
    splitt=np.zeros((numrs*2,boxesperR+2))
    while(index<len(H)):
        if(float(H[index])!=-3.0*3.96):
            #print(H[index])
            array=np.append(array,[float(H[index])],axis=0)
            index+=1
        else:
            print "Gonna make a histogram"
            array=array[1:]
            v=u.histogram(array,int(boxesperR))
            v=v[1:]
            Ni=0
            for c in range(0,len(v)):
                Ni+=v[c][1]
            print Ni
            for g in range(0,len(v)):
                splitt[2*rnum][g]=v[g][0]
                splitt[2*rnum+1][g]=v[g][1]
            rnum+=1
            print(len(array))
            array=np.array([0.0])
            index+=1
         
    return splitt

###used to get the histogram for each window with a certain box width
def Hsplitwidth(Hfile,numrs,width):
    with open(Hfile) as f:
        data = f.readlines()
        #data= data[199328:]
    H=np.array([-3.0])
    for row in data:
        H=np.append(H,[float(row.split()[0])*3.96],axis=0)
    index=0
    aindex=0
    rnum=0
    H=H[1:]
    #print 'H array made'
    min=10000
    max=0
    widthmax=0
    notsplitt=np.zeros((numrs,(len(H)/numrs)-1))
    while(index<len(H)):
        if(float(H[index])!=-3.0*3.96):
            #print(H[index])
            #print rnum
            #print aindex
            notsplitt[rnum][aindex]=float(H[index])
            if(float(H[index])>max):
                max=float(H[index])
            if float(H[index]<min):
                min=float(H[index])
            index+=1
            aindex+=1
        else:
            if((max-min)>widthmax):
                widthmax=max-min
                print max
                print min
            rnum+=1
            index+=1
            aindex=0
            min=10000
            max=0
    print width
    print widthmax    
    splitt=np.zeros((numrs*2,(int(widthmax/width)+1)))
    print "Gonna make a histogram"
    for z in range(0,numrs):
        max=0
        min=10000
        for y in range(0,len(notsplitt[z])):
            if(notsplitt[z][y]>max):
                max=notsplitt[z][y]
            if (notsplitt[z][y]<min):
                min=notsplitt[z][y]
                print min
        boxes=int((max-min)/width)-2
        v=u.histogram(notsplitt[z],int(boxes))
        v=v[1:]
        print boxes
        print splitt.shape
        print v.shape
        for c in range(0,boxes+2):
            splitt[2*z][c]=v[c][0]
            print  splitt[2*z][c]
            splitt[2*z+1][c]=v[c][1]
            print splitt[2*z+1][c]
    return splitt
        
def H_not_split(Hfile,numrs,boxesperR):
    with open(Hfile) as f:
        data = f.readlines()
        data=data[9309:]
    H=np.array([-3.0])
    
    for row in data:
        if(float(row.split()[0])!=-3.0):
            H=np.append(H,[float(row.split()[0])*3.96],axis=0)
    H=H[1:] 
    print H[1]
    v=u.histogram(H,int(boxesperR*numrs))
    v=v[1:]
    return v


#makes one total histogram with a ceratin box width
def H_not_split_width(Hfile,numrs,width):
    with open(Hfile) as f:
        data = f.readlines()
       # data=data[9309:]
    H=np.array([-3.0])
    max=0
    min=1000
    for row in data:
        if(float(row.split()[0])!=-3.0):
            H=np.append(H,[float(row.split()[0])*3.96],axis=0)
            if(float(row.split()[0])<min):
                min=float(row.split()[0])
                print min
            if(float(row.split()[0])>max):
                max=float(row.split()[0])
                print max
    H=H[1:]
    print ""
    print max
    print min
    print max-min
    print (max-min)*3.96/width
    print int((max-min)/width)
    v=u.histogram(H,int((max-min)*3.96/width))
    v=v[1:]
    return v


