# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""

import sys
sys.path.insert(0,'../analysis/hydrocarbon')
import numpy as np
import matplotlib.pyplot as plt
import util as u
#import pair_force_solver as pfs
import pfsupload as pfs

boxesperR=20
windows=11
width=.1
array=pfs.Hsplitwidth('H.txt',windows,width)
u.write_array(array,'output.txt')
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
Ni=1
#8.18+1.97+1.25
npRadius=1.14
repeats=12
l=.128*repeats+.2
lamda=l/npRadius
opmtau=(3*lamda+1)**(1/3)
ocmtau=-1*(1+lamda)/2+(((1+lamda)/2)**2+(6*lamda+2)/(1+lamda))**(.5)


with open("output.txt") as f:
    data = f.readlines()
print data[0].split()[0]
allR=np.array([0])
for i in range(0,len(data)/2):
    for x in range(0,len(data[0].split())):
        allR=np.append(allR,[float(data[2*i].split()[x])],axis=0)
windows=len(data)/2
for w in range(0,windows):
    y=np.zeros(len(allR))
    for g in range(0,len(data[1].split())):
        for x in range(0,len(allR)):
            if (float(data[2*w].split()[g])==allR[x]):
                y[x]=float(data[2*w+1].split()[g])/Ni
    print y
    ax1.plot(allR,y, c=[(w%3==0),(w%3==2),(w%3==1)], label=w, linewidth=1)

spacing=.5
for t in range(0,windows):
	ax1.plot([75-t*2],[0],marker='.',c=[(t%3==0),(t%3==2),(t%3==1)],markersize=20,label=str(t)+'th bias')              
""""  
ax1.plot([57.5,57.5,[0,5000],c='0.0',label='first interactions')
ax1.plot([45,45],[0,5000],c='0.4',label='first touch')
ax1.plot([20,20],[0,5000],c='0.8',label='sulfur touch')
ax1.plot([47.5,47.5],[0,5000],c='0.2',label='first inflection')
ax1.plot([29.5,29.5],[0,5000],c='0.6',label='seond inflection')
"""
ax1.plot([2*npRadius*10*opmtau],[0], linestyle='None', marker='.', color='red', markersize=20,label='OPMcut')
ax1.plot([2*npRadius*10*ocmtau],[0], linestyle='None', marker='.', color='orange', markersize=20,label='OCMcut')
ax1.set_title('Denstity')    
ax1.set_xlabel('R')
ax1.set_ylabel('Density')
    
leg1 = ax1.legend()

plt.show()

