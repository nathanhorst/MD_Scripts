"""
Created on Mon Jul 11 15:04:51 2016

@author: waltmann
"""


import numpy as np

"""
takes a list and returns the normla distribution of the list
not currently useful, but who knows.
"""
def normal_distribution(arr):
    std=np.std(arr)
    mean=np.mean(arr)
    dist=np.array([[0.0,0.0]])
    d=1/np.sqrt(2*std*std*np.pi)
    for i in range(0,60): 
        x=(mean-3*std)+i*(std/10.0)
        a=np.exp(-1*((x-mean)**2)/(2*std**2))*d
        dist=np.append(dist,[[x,a]],axis=0)
    return dist    
    
#print(len(inter_angle('2np.xml',0)))  


"""
creates a histgoram of some list given a number of boxes. make sure boxes is  an integer
"""    
def histogram(arr,boxes):    
    arr=np.sort(arr,axis=0,kind='mergesort')
    #print arr
    hist=np.array([[0.0,0.0]])
    step= (arr[len(arr)-1]-arr[0])/boxes
    for i in range(0,boxes+2):
        min=arr[0]+(i-1)*step
        max=arr[0]+(i)*step
        hit=0
        for x in range(0,len(arr)):
            if (arr[x]>=min and arr[x]<max):
                hit+=1
        hist=np.append(hist,[[(min+max)/2,hit]],axis=0)
    #print hist
    return hist
    
def vector_length(x,y,z):
    c=np.sqrt(x**2+y**2+z**2)
    return c
 

"""
tkaes a vetor and returns its length
"""   
def length(vector):
    return vector_length(vector[0],vector[1],vector[2])   
    
"""
calculates the distance between two particles even if one crossed the box.
Must know box dimensions. could be more useful than the file below if box dimensions
are known and you don't want to keep reading the same file
"""
def part_distance_box_dimensions(part1,part2,lx,ly,lz):
    
    x=abs(float(part1[0])-float(part2[0]))
    if(x>(lx/2)-1):
        x=lx-x
    y=abs(float(part1[1])-float(part2[1]))
    if(y>(ly/2)-1):
        y=ly-y
    z=abs(float(part1[2])-float(part2[2]))
    if(z>(lz/2)-1):
        z=lz-z
    return vector_length(x,y,z)

"""
does the distance between particles without knowing the box size expliciitly.
just the file it came from
"""
def part_distance(part1,part2,inputfile):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    #print data1[3]
    s=data1[3].split()
    #print s
    lx=0
    ly=0
    lz=0
    xy=0
    for i in range(len(s)):
        if s[i][:2]=='lx':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            lx=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='ly':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            ly=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='lz':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            lz=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='xy':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            xy=float(s[i][count+1:b+count-len(s[i])])
    #print lx
    #print ly
    #print lz
    #print part1
    #print part2[0]
    x=abs(float(part1[0])-float(part2[0]))
    if(x>(lx/2)-1):
        x=lx-x
    y=abs(float(part1[1])-float(part2[1]))
    
    z=abs(float(part1[2])-float(part2[2]))
    if(z>(lz/2)-1):
        z=lz-z
    #print x
    #print y 
    #print z
    vector=[x,y,z]
    #print vector
    v2=[xy*ly,ly,0]
    v2hat=v2/length(v2)
    ycheck=np.dot((vector),v2hat)
    #print('v2hat is '+ str(v2hat))
    #print np.subtract(part1,part2)
    #print('ycheck is '+ str(ycheck))
    if(ycheck>(length(v2)/2)-1):
        vector=np.subtract(vector,v2)
        print vector
    if(vector[0]>(lx/2)-1):
        vector[0]=lx-vector[0]
        print vector
    return vector_length(vector[0],vector[1],vector[2])
#print part_distance([-70,-39,39],[-39,39,-39],'Lat108_25_Fcc.xml') 

"""
instead of the distance it returns the vector from part1 to part2
"""    
def part_vector(part1,part2,inputfile):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    #print data1[3]
    s=data1[3].split()
    #print s
    lx=0
    ly=0
    lz=0
    for i in range(len(s)):
        if s[i][:2]=='lx':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            lx=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='ly':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            ly=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='lz':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            lz=float(s[i][count+1:b+count-len(s[i])])
    #print lx
    #print ly
    #print lz
    #print part1
    #print part2[0]
    x=abs(float(part1[0])-float(part2[0]))
    if(x>(lx/2)-1):
        x=lx-x
    y=abs(float(part1[1])-float(part2[1]))
    #print y
    if(y>(ly/2)-1):
        y=ly-y
    z=abs(float(part1[2])-float(part2[2]))
    if(z>(lz/2)-1):
        z=lz-z
    #print x
    #print y 
    #print z
    l=np.array([x,y,z])
    return l
   

"""
used by many things to write to text files.
always cuts off the first row (index 0) of what is assumed to be a 2d array
"""
def write_array(arr, outputfile):
    fid=open(outputfile,'w')
    for i in range(1,arr.shape[0]):
        for x in range(0,arr.shape[1]):
            fid.write(str(arr[i][x]))
            if(x!=arr.shape[1]-1):
                fid.write(" ")
            else:
                fid.write('\n')


"""
used by radial distribution
changing a in the first line can change the sotness of final curves. as a approaches
0 it becomes a truer dirac delta funciton
"""
def dirac_delta(r,lower_limit,upper_limit,step):
    a=.1
    b=(upper_limit-lower_limit)/float(step)
    result=np.array([[0.0,0.0]])
    for i in range(0,step):
        n=lower_limit+b*i
        #print n
        y=(1/a*np.sqrt(np.pi))*np.exp(-1*((n-r)**2/a**2))
        result=np.append(result,[[n,y]],axis=0)
    #print result
    return result
#write_array(dirac_delta(1.9343,0,2.4234234,200),'dd.txt')     
    
    
"""
returns an array of the box dimensions
"""    
def read_box_dimensions(inputfile):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    #print data1[3]
    s=data1[3].split()
    #print s
    lx=0
    ly=0
    lz=0
    for i in range(len(s)):
        if s[i][:2]=='lx':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            lx=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='ly':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            ly=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='lz':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            lz=float(s[i][count+1:b+count-len(s[i])])
    return np.array([lx,ly,lz])
    
