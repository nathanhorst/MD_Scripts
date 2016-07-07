# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
import numpy as np
import matplotlib.pyplot as plt

with open("trans_Tvar.txt") as f:
    data = f.readlines()

with open("trans_Tvar1.txt") as f1:
    data1 = f1.readlines()
with open("trans_Tvar2.txt") as f2:
    data2 = f2.readlines()
with open("trans_Tvar3.txt") as f3:
    data3 = f3.readlines()
with open("trans_Tvar4.txt") as f4:
    data4 = f4.readlines()
with open("trans_Tvar5.txt") as f5:
    data5 = f5.readlines()
with open("trans_Tvar6.txt") as f6:
    data6 = f6.readlines()
with open("trans_Tvar7.txt") as f7:
    data7 = f7.readlines()
with open("trans_Tvar8.txt") as f8:
    data8 = f8.readlines()
with open("trans_Tvar9.txt") as f9:
    data9 = f9.readlines()
with open("trans_Tvar10.txt") as f10:
    data10 = f10.readlines()

x = [float(row.split()[0]) for row in data]
y1 = [float(row.split()[1])*100.0 for row in data]
y2 = [float(row.split()[1])*100.0 for row in data1]
y3 = [float(row.split()[1])*100.0 for row in data2]
y4 = [float(row.split()[1])*100.0 for row in data3]
y5 = [float(row.split()[1])*100.0 for row in data4]
y6 = [float(row.split()[1])*100.0 for row in data5]
y7 = [float(row.split()[1])*100.0 for row in data6]
y8 = [float(row.split()[1])*100.0 for row in data7]
y9 = [float(row.split()[1])*100.0 for row in data8]
y10 = [float(row.split()[1])*100.0 for row in data9]
y11 = [float(row.split()[1])*100.0 for row in data10]

fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

var1='Percent Trans (20-500K)'

ax1.set_title(var1)    
ax1.set_xlabel('Timestep')
ax1.set_ylabel(var1)

ax1.plot(x,y1, c='black', label='all',linewidth=3)
ax1.plot(x,y2, c='red', label='1st',linewidth=1)
ax1.plot(x,y3, c='orange', label='2nd',linewidth=1)
ax1.plot(x,y4, c='yellow', label='3rd',linewidth=1)
ax1.plot(x,y5, c='green', label='4th',linewidth=1)
ax1.plot(x,y6, c='blue', label='5th',linewidth=1)
ax1.plot(x,y7, c='cyan', label='6th',linewidth=1)
ax1.plot(x,y8, c='purple', label='7th',linewidth=1)
ax1.plot(x,y9, c='brown', label='8th',linewidth=1)
ax1.plot(x,y10, c='gray', label='9th',linewidth=1)
ax1.plot(x,y11, c='pink', label='10th',linewidth=1)


leg1 = ax1.legend()

plt.show()
