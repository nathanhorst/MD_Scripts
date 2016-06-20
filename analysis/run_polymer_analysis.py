# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 09:04:02 2016

@author: waltmann
"""

import os
import sys

import numpy as np
from numpy import linalg as LA

import torsional as t
import distance_nth_carbon as dnc


########
inputfile='atoms2.dump.0004500000.xml'

total=0.0

for i in range(1,13):
    q=dnc.distance_to_nth_carbon(inputfile,i)
    for f in range(len(q)):
        total+=q[f]
    print (str(i)+'  '+str(total/(float(f)+1)))
    total=0

numtrans=0.0
numgplus=0.0
numgminus=0.0    
trans_cut=1.0472
d=t.all_dihedral_angles(inputfile)
for x in range(1,len(d)):
    if(abs(d[x])<trans_cut):
        numtrans+=1;
    elif(d[x]<-1.0*trans_cut):
        numgminus+=1
    else:
        numgplus+=1
print('The percent trans is ' + str(numtrans/float(len(d)-1)))
print('The percent gauche minus is ' + str(numgminus/float(len(d)-1)))
print('The percent gauche plus is ' + str(numgplus/float(len(d)-1)))

