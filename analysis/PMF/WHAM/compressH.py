# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 11:05:53 2017

@author: waltmann
"""

import numpy as np

import pfsupload as pfs

def even_H(files):
    min=10e45
    print min
    print len(files)
    for i in range(0,len(files)):
        fcn=open(files[i],'r')
	print files[i]
	x=fcn.readlines()
        if (len(x)<min):
	    print len(x)
            min=len(x)
	    print min 	
    print min
    for i in range(0,len(files)):
        fcn=open(files[i],'r')
        pfs.cut_lines_H_end(files[i],len(fcn.readlines())-min,'new'+files[i])
	pfs.cut_lines_H('new'+files[i],1,files[i])
#['H_123.txt','H_126.txt']
def compress_h(hs):
    #strin=np.chararray(len(hs),itemsize=15)    
    num=np.zeros(len(hs))
    for y in range(len(hs)):
        print float(hs[y][2:-4])
        #strin[y]=hs[y][2:-4]
        num[y]=float(hs[y][2:-4])
    #num = num[::-1]
    #print num
    #num=np.sort(num,kind='mergesort')
    r=open('R.txt','w+')
    h=open('H.txt','w+')
    for i in range(len(hs)):
        r.write(str(num[i]))
        r.write('\n')
        temp=open('H_'+ str(hs[i][2:-4])+ '.txt','r')
        data=temp.readlines()
        for line in data:
            h.write(line)
        h.write('-3.0')
        h.write('\n')
        
t=[25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59]
files=len(t)
b=np.chararray(files,itemsize=15)
for i in range(files):
    v='H_'+ str(t[files-i-1]) + '.0.txt'
    print v
    b[i]=v
print b
even_H(b)
compress_h(b)
