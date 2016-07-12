# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
import numpy as np
import matplotlib.pyplot as plt

with open("mylog.log") as f:
    data = f.readlines()
    data = data[2:]

x = [row.split()[0] for row in data]
y1 = [row.split()[1] for row in data]
y2 =[row.split()[2] for row in data]
y3 = [float(row.split()[3])*(2/90030.) for row in data]
y4 =[row.split()[4] for row in data]
y5 =[row.split()[5] for row in data]

fig1 = plt.figure()
#fig2 = plt.figure()
#fig3 = plt.figure()
#fig4 = plt.figure()
#fig5 = plt.figure()

ax1 = fig1.add_subplot(111)
#ax2 = fig2.add_subplot(111)
#ax3 = fig3.add_subplot(111)
#ax4 = fig4.add_subplot(111)
#ax5 = fig5.add_subplot(111)

var1='Temp and KE'
#var2='Potential Energy'
#var3='Kinetic Energy'
#var4='Volume'
#var5='Pressure'

ax1.set_title(var1)    
ax1.set_xlabel('Timestep')
ax1.set_ylabel(var1)
#ax2.set_title(var2)    
#ax2.set_xlabel('Timestep')
#ax2.set_ylabel(var2)
#ax3.set_title(var3)    
#ax3.set_xlabel('Timestep')
#ax3.set_ylabel(var3)
"""ax4.set_title(var4)    
ax4.set_xlabel('Timestep')
ax4.set_ylabel(var4)
ax5.set_title(var5)    
ax5.set_xlabel('Timestep')
ax5.set_ylabel(var5)
"""
ax1.plot(x,y1, c='r', label='T')
#ax2.plot(x,y2, c='blue', label='PE')
ax1.plot(x,y3, c='g', label='KE')
#ax4.plot(x,y4, c='b', label='V')
#ax5.plot(x,y5, c='orange', label='P')

leg1 = ax1.legend()
#leg2 = ax2.legend()
#leg3 = ax3.legend()
#leg4 = ax4.legend()
#leg5 = ax5.legend()

plt.show()
