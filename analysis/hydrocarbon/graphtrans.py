# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:32:15 2016

@author: nathan
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors

def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv') 
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

length=10
timesteps=100
data=np.zeros(((2*length-1),timesteps))
for i in range(0,length-2):
	with open("trans_Tvar"+ str(i) +".txt") as f1:
    		data1 = f1.readlines()
		data[2*(i)]=[float(row.split()[0]) for row in data1]
		data[2*(i)+1]= [float(row.split()[1])for row in data1]


for i in range(0,length-2):
    total=0
    runs=0
    for g in range(0,timesteps):
        total+=data[2*i+1][g] 
    print 'average trans ' + str(i) + ' is ' +str((total/timesteps))






fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

var1='Percent Trans'

ax1.set_title(var1)    
ax1.set_xlabel('Timestep')
ax1.set_ylabel(var1)

for d in range(0,length-2):
	x=1
	if(d==0):
		x=3
	ax1.plot(data[2*d],data[2*d+1], c=np.random.rand(3,), label=str(d),linewidth=x)


leg1 = ax1.legend()

plt.show()
