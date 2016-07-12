# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
import numpy as np
import matplotlib.pyplot as plt

with open("density.txt") as f:
    data = f.readlines()

width=0.25

x = [float(row.split()[0])*3.96 for row in data]
y1 = [float(row.split()[1])/(3.96**3) for row in data]

fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

var1='Particle Density'

ax1.set_title(var1)    
ax1.set_xlabel('Distance (Ang)')
ax1.set_ylabel(var1)

ax1.bar(x,y1,width, color='black')

leg1 = ax1.legend()

plt.show()
