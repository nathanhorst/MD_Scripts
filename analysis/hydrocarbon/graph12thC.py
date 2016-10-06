# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
import numpy as np
import matplotlib.pyplot as plt


length=19
timesteps=100
data=np.zeros(((2*length-1),timesteps))
for i in range(1,length):
	with open("dist"+ str(i) +"thC.txt") as f1:
    		data1 = f1.readlines()
		data[2*(i-1)]=[float(row.split()[0]) for row in data1]
		data[2*(i-1)+1]= [float(row.split()[1])*3.96 for row in data1]


for i in range(0,length-1):
    total=0
    runs=0
    for g in range(0,timesteps):
        total+=data[2*i+1][g] 
    print 'average distance' + str(i+1) + ' is ' +str((total/timesteps))






fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

var1='Distance to nth C (Ang)'

ax1.set_title(var1)    
ax1.set_xlabel('Timestep')
ax1.set_ylabel(var1)

for d in range(0,length-1):
	ax1.plot(data[2*d],data[2*d+1], c=[.05+d*(1.0/length),0,0], label=str(d),linewidth=1)


leg1 = ax1.legend()

plt.show()
