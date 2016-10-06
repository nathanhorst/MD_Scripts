# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
import numpy as np
import matplotlib.pyplot as plt

with open("dist_0_1.txt") as f1:
    data1 = f1.readlines()
with open("dist_1_2.txt") as f2:
    data2 = f2.readlines()
with open("dist_2_3.txt") as f3:
    data3 = f3.readlines()
with open("dist_3_4.txt") as f4:
    data4 = f4.readlines()
with open("dist_4_5.txt") as f5:
    data5 = f5.readlines()
with open("dist_5_6.txt") as f6:
    data6 = f6.readlines()
with open("dist_6_7.txt") as f7:
    data7 = f7.readlines()
with open("dist_7_8.txt") as f8:
    data8 = f8.readlines()
with open("dist_8_9.txt") as f9:
    data9 = f9.readlines()
with open("dist_9_10.txt") as f10:
    data10 = f10.readlines()
with open("dist_10_11.txt") as f11:
    data11 = f11.readlines()
with open("dist_11_12.txt") as f12:
    data12 = f12.readlines()

x = [float(row.split()[0]) for row in data1]
y1 = [float(row.split()[1])*3.96 for row in data1]
y2= [float(row.split()[1])*3.96 for row in data2]
y3= [float(row.split()[1])*3.96 for row in data3]
y4= [float(row.split()[1])*3.96 for row in data4]
y5= [float(row.split()[1])*3.96 for row in data5]
y6= [float(row.split()[1])*3.96 for row in data6]
y7= [float(row.split()[1])*3.96 for row in data7]
y8= [float(row.split()[1])*3.96 for row in data8]
y9= [float(row.split()[1])*3.96 for row in data9]
y10= [float(row.split()[1])*3.96 for row in data10]
y11= [float(row.split()[1])*3.96 for row in data11]
y12= [float(row.split()[1])*3.96 for row in data12]



q=np.array([data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12])

for i in range(0,len(q)):
    total=0
    runs=0
    for g in range(0,len(data1)):
        total+=float(q[i][g].split()[1])
        runs+=1 
    print 'average distance' + str(i+1) + ' is ' +str((total/runs)*3.96)






fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

var1='Distance to nth C (Ang)'

ax1.set_title(var1)    
ax1.set_xlabel('Timestep')
ax1.set_ylabel(var1)

ax1.plot(x,y1, c='black', label='1',linewidth=1)
ax1.plot(x,y2, c='red', label='2',linewidth=1)
ax1.plot(x,y3, c='orange', label='3',linewidth=1)
ax1.plot(x,y4, c='yellow', label='4',linewidth=1)
ax1.plot(x,y5, c='green', label='5',linewidth=1)
ax1.plot(x,y6, c='blue', label='6',linewidth=1)
ax1.plot(x,y7, c='purple', label='7',linewidth=1)
ax1.plot(x,y8, c='pink', label='8',linewidth=1)
ax1.plot(x,y9, c='brown', label='9',linewidth=1)
ax1.plot(x,y10, c='gray', label='10',linewidth=1)
ax1.plot(x,y11, c='maroon', label='11',linewidth=1)
ax1.plot(x,y12, c='cyan', label='12',linewidth=1)

leg1 = ax1.legend()

plt.show()
