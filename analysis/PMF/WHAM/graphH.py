import numpy as np
import matplotlib.pyplot as plt
import correlation as cor













#r=cor.estimated_autocorrelation(cor.get_positions('18.txt'))
r=cor.get_positions('H_69.0.txt')

print np.mean(r[:2500])

print np.mean(r[2500:])
fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

x=np.linspace(1,len(r),len(r))

var1='Correlation'
var2='K'
#var3='Kinetic Energy'
#var4='Volume'
#var5='Pressure'

ax1.set_title(var1)    
ax1.set_xlabel('K') 



ax1.plot(x,r, c='r', label='Corr')

leg1 = ax1.legend()

plt.show()
