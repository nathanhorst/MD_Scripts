# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
import numpy as np
import matplotlib.pyplot as plt

with open("dd.txt") as f:
    data = f.readlines()
    data = data[1:]
peak=0
length=0
for row in data:
	s=row.split()
	if (float(s[1])>peak):
		peak=float(s[1])
		length=float(s[0])

#print length
#print peak

x = [(float(row.split()[0])) for row in data]
y1 = [(float(row.split()[1]))/peak for row in data]

fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

var1='Radial Density'

ax1.set_title(var1)    
ax1.set_xlabel('spacing')
ax1.set_ylabel(var1)
#ax3.set_title(var3)    
#ax3.set_xlabel('Timestep')
#ax3.set_ylabel(var3)
ax1.plot(x,y1, c='r', label='Intensity')

leg1 = ax1.legend()

plt.show()
