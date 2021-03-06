# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

with open("mylog.log") as f:
    data = f.readlines()
    data = data[2:]

npart=1
ncore=123
nchains=80
L=22
rigid=npart
nonrigid=nchains*(L-1)

N=npart*(rigid+nonrigid)

x = [row.split()[0] for row in data]
y1 = [row.split()[1] for row in data]
y2 =[float(row.split()[2]) for row in data]
y3 = [float(row.split()[3])*(2./(3.*N)) for row in data]
y4 =[row.split()[4] for row in data]
y5 =[row.split()[5] for row in data]
y6 =[float(row.split()[6]) for row in data]
y7 =[float(row.split()[7]) for row in data]

fig1 = plt.figure()
fig2 = plt.figure()
#fig3 = plt.figure()
#fig4 = plt.figure()

ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
#ax3 = fig3.add_subplot(111)
#ax4 = fig4.add_subplot(111)


var1='Temp and KE'
var2='Energies'
#var3='Volume'
#var4='Pressure'

ax1.set_title(var1)    
ax1.set_xlabel('Timestep')
ax1.set_ylabel(var1)
ax2.set_title(var2)    
ax2.set_xlabel('Timestep')
ax2.set_ylabel(var2)
#ax3.set_title(var3)    
#ax3.set_xlabel('Timestep')
#ax3.set_ylabel(var3)
#ax4.set_title(var4)    
#ax4.set_xlabel('Timestep')
#ax4.set_ylabel(var4)

ax1.plot(x,y1, c='r', label='T')
ax2.plot(x,y2, c='black', label='PE')
ax1.plot(x,y3, c='g', label='KE')
#ax3.plot(x,y4, c='b', label='V')
#ax4.plot(x,y5, c='orange', label='P')
#ax2.plot(x,y6, c='r', label='FS_LJ')
#ax2.plot(x,y7, c='blue', label='bond_fene')

leg1 = ax1.legend()
leg2 = ax2.legend()
#leg3 = ax3.legend()
#leg4 = ax4.legend()

plt.show()
