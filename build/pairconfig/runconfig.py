# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:52:04 2016

@author: nathan
"""
from __future__ import division
        
def build_2NP_PF_xyz(read,save,R):

    f2=open(read,'r')
    f=open(save,'w')
    fdata=f2.read()
    fd=fdata.splitlines()
    N=fd[0]
    print fd[0]  
    fd=fd[2:]
    print fd[0]
 
    f.write(('%s\n\n')%(str(2*int(N))))
    for line in fd:
        s=line.split()
        x=float(s[1])-float(R/2)
        y=float(s[2])
        z=float(s[3])
	if s[0]=='V1':
		f.write(("V1 %f %f %f\n")%(x,y,z))
	else:
        	f.write(("%s %f %f %f\n")%(s[0],x,y,z))
    for line in fd:
        s=line.split()
        x=float(s[1])+float(R/2)
        y=float(s[2])
        z=float(s[3])
	if s[0]=='V1':
        	f.write(("V2 %f %f %f\n")%(x,y,z))
	else:
		f.write(("%s %f %f %f\n")%(s[0],x,y,z))

build_2NP_PF_xyz('PEOsphere29.xyz','PEOsphere29_r69.xyz',69)
