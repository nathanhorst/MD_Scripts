
import numpy as np
import matplotlib.pyplot as plt

with open("pressure.log") as f:
    data = f.readlines()
    data = data[2:]

x = [row.split()[0] for row in data]
y1 = [row.split()[1] for row in data]
y2 =[row.split()[2] for row in data]
y3 = [row.split()[3] for row in data]
y4 =[row.split()[4] for row in data]
y5 =[row.split()[5] for row in data]
y6 =[row.split()[6] for row in data]
y7 =[row.split()[7] for row in data]


fig = plt.figure()

ax1 = fig.add_subplot(111)

var='Pressure'

ax1.set_title(var)    
ax1.set_xlabel('Timestep')
ax1.set_ylabel(var)

ax1.plot(x,y1, c='black',linewidth=2, label='P',)
ax1.plot(x,y2, c='blue',linewidth=0.5, label='Pxx')
ax1.plot(x,y3, c='green',linewidth=0.5, label='Pyy')
ax1.plot(x,y4, c='red',linewidth=0.5, label='Pzz')
ax1.plot(x,y5, c='orange',linewidth=0.5, label='Pxy')
ax1.plot(x,y6, c='purple',linewidth=0.5, label='Pyz')
ax1.plot(x,y7, c='yellow',linewidth=0.5, label='Pxz')


leg = ax1.legend(bbox_to_anchor=(1.1,1.05))

plt.show()
