# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
import numpy as np
import matplotlib.pyplot as plt

with open("mylog.log") as f:
    data = f.readlines()
    data = data[20:]

x = [row.split()[4] for row in data]
y1 = [(float(row.split()[5])) for row in data]

fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

var1='Pressure'

ax1.set_title(var1)    
ax1.set_xlabel('Volume')
ax1.set_ylabel(var1)
#ax3.set_title(var3)    
#ax3.set_xlabel('Timestep')
#ax3.set_ylabel(var3)
ax1.plot(x,y1, c='r', label='Pressure')

leg1 = ax1.legend()

plt.show()
