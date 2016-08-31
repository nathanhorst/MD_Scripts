# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
import numpy as np
import matplotlib.pyplot as plt

with open("dist_end.txt") as f:
    data = f.readlines()

x = [float(row.split()[0]) for row in data]
y1 = [float(row.split()[1]) for row in data]

fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

var1='End to End Distance (Ang)'

ax1.set_title(var1)    
ax1.set_xlabel('Timestep')
ax1.set_ylabel(var1)

ax1.plot(x,y1, c='black', label='all',linewidth=3)

leg1 = ax1.legend()

plt.show()
