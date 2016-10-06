###################################################################
## This file can be used to graph the output log file from HOOMD ##
###################################################################
import numpy as np
import matplotlib.pyplot as plt

# This statement opens the file and read in the data, omitting the first two lines of the simulation
with open("mylog.log") as f:
    data = f.readlines()
    data = data[2:]

# The following splits up the data read from the log file by column, setting the values in the first column equal to the x coordinate and the values in the second column to the y coordinate. If you output more quantities in the log files, you can set additional variables y2, y3, etc. to read the corresponding columns by changing the number in the [] accordingly
x = [row.split()[0] for row in data]
y1 = [row.split()[1] for row in data]

# The following lines set the axis labels, legend location, etc. and will need to be modified for each additional y variable you impose. Take a look at the matplotlib documentation if you are unsure what these variables are indicating.
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
var1='Potential Energy'
ax1.set_title(var1)    
ax1.set_xlabel('Timestep')
ax1.set_ylabel(var1)
ax1.plot(x,y1, c='k', label='PE')
leg1 = ax1.legend()
plt.show()
