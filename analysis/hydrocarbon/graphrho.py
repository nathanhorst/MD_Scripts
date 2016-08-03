# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
import numpy as np
import matplotlib.pyplot as plt
import radial_distribution as rad
import util as u



with open("H.txt") as f:
    data = f.readlines()
H=np.array([-3.0])
for row in data:
    H=np.append(H,[float(row.split()[0])],axis=0)
H=H[1:]
print H
v=u.histogram(H,int(10*3.96*10))
v=v[1:]
print v
x = [float(row[0])*3.96 for row in v]
y1 = [float(row[1])/len(H) for row in v]


fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

var1='probability'

ax1.set_title(var1)    
ax1.set_xlabel('r')
ax1.set_ylabel(var1)

ax1.plot(x,y1, c='black',label='all',linewidth=1)
width=.0001
#ax1.bar(x,y1,width, color='black')

plt.show()
