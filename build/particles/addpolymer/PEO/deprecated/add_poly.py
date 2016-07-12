# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import sys
import os
from numpy import *
import numpy as np

#def poly_from_points(inputfile,outputfile,R1,R2,R3,L):
def poly_from_points(inputfile,outputfile,R1,R2,L):    
    #######################################################################################################################################################
    ### This function creates an xyz file including polymers of length L extending off points of the crystal at points that are distance R from the origin
    #######################################################################################################################################################
    
    #Open xyz input file to append with new points
    
    fin=open(inputfile,'r')
    fout=open(outputfile,'w')
    fdata=fin.read()
    data=fdata.splitlines()
    data=data[2:]    
    
    for line in data:
        s=line.split()
        fout.write(('%s %f %f %f\n')%(s[0],float(s[1]),float(s[2]),float(s[3])))
        
    j=0    
    for line in data:
        ## shift value determines spacing of polymer addition
        shift=0.3
        s=line.split()
        vec=np.sqrt((float(s[1])**2)+(float(s[2])**2)+(float(s[3])**2))        
        if vec == R1:
            j += 1
            print(j)
            for i in range(L):
                if i == 0:                    
                    r1=float(s[1])+shift*float(s[1])
                    r2=float(s[2])+shift*float(s[2])
                    r3=float(s[3])+shift*float(s[3])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(s[1])+shift*(i+1)*float(s[1])
                    r2=float(s[2])+shift*(i+1)*float(s[2])
                    r3=float(s[3])+shift*(i+1)*float(s[3])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
        if vec == R2:
            j += 1
            print(j)
            for i in range(L):
                if i == 0:                    
                    r1=float(s[1])+shift*float(s[1])
                    r2=float(s[2])+shift*float(s[2])
                    r3=float(s[3])+shift*float(s[3])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(s[1])+shift*(i+1)*float(s[1])
                    r2=float(s[2])+shift*(i+1)*float(s[2])
                    r3=float(s[3])+shift*(i+1)*float(s[3])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))

'''                    
Example script to run simulation -- builds fccnp_cutbcc with polymers on the 111 faces 

import sys
import os
import numpy as np

import add_poly as ap

#################
## File to read
#################

inputfile='bccnp.xyz'
outputfile='bccnp_withpoly.xyz'

#############################################
# Define coordinates of points to build from
#############################################

x1=4
y1=1
z1=0

R1=np.sqrt(x1**2+y1**2+z1**2)

###########

x2=2
y2=2
z2=2

R2=np.sqrt(x2**2+y2**2+z2**2)

###################################
# Length of polymer chain to build
###################################

L=10

##############
# Run Builder
##############

ap.poly_from_points(inputfile,outputfile,R1,R2,L)
                 
'''

                    
def poly_without_points(inputfile,outputfile,R1,L,shift1):

    #######################################################################################################################################
    ### This function appends the xyz file made from poly_from_points to include polymers of length L on the 001 faces of the fccnp_cutbcc
    #######################################################################################################################################   
    
    ###############################################################################################################################
    ### !!!!!!!!!CAUTION!!!!!!!!!!!! where the previous function is general to build polymers from points of a specified length
    ### from the center of the particle, this function is highly geometry dependent and needs adjusting for different shaped particles    
    ###############################################################################################################################
    
    
    fin=open(inputfile,'r')
    fout=open(outputfile,'w')
    fdata=fin.read()
    data=fdata.splitlines()
    data=data[2:]    
    
    for line in data:
        s=line.split()
        fout.write(('%s %f %f %f\n')%(s[0],float(s[1]),float(s[2]),float(s[3])))
    
    for line in data:
        shift=0.3
        s=line.split()
        vec=np.sqrt((float(s[1])**2)+(float(s[2])**2)+(float(s[3])**2))
        if vec == R1 and float(s[1]) != float(0):
            pos1=[float(s[1]),float(s[2])+float(shift1),float(s[3])]
            for i in range(L):
                if i == 0:
                    r1=float(pos1[0])+shift*float(pos1[0])
                    r2=float(pos1[1])+shift*float(pos1[1])
                    r3=float(pos1[2])+shift*float(pos1[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos1[0])+shift*(i+1)*float(pos1[0])
                    r2=float(pos1[1])+shift*(i+1)*float(pos1[1])
                    r3=float(pos1[2])+shift*(i+1)*float(pos1[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
            pos2=[float(s[1]),float(s[2])-float(shift1),float(s[3])]
            for i in range(L):
                if i == 0:
                    r1=float(pos2[0])+shift*float(pos2[0])
                    r2=float(pos2[1])+shift*float(pos2[1])
                    r3=float(pos2[2])+shift*float(pos2[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos2[0])+shift*(i+1)*float(pos2[0])
                    r2=float(pos2[1])+shift*(i+1)*float(pos2[1])
                    r3=float(pos2[2])+shift*(i+1)*float(pos2[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
            pos3=[float(s[1]),float(s[2]),float(s[3])+float(shift1)]
            for i in range(L):
                if i == 0:
                    r1=float(pos3[0])+shift*float(pos3[0])
                    r2=float(pos3[1])+shift*float(pos3[1])
                    r3=float(pos3[2])+shift*float(pos3[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos3[0])+shift*(i+1)*float(pos3[0])
                    r2=float(pos3[1])+shift*(i+1)*float(pos3[1])
                    r3=float(pos3[2])+shift*(i+1)*float(pos3[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
            pos4=[float(s[1]),float(s[2]),float(s[3])-float(shift1)]
            for i in range(L):
                if i == 0:
                    r1=float(pos4[0])+shift*float(pos4[0])
                    r2=float(pos4[1])+shift*float(pos4[1])
                    r3=float(pos4[2])+shift*float(pos4[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos4[0])+shift*(i+1)*float(pos4[0])
                    r2=float(pos4[1])+shift*(i+1)*float(pos4[1])
                    r3=float(pos4[2])+shift*(i+1)*float(pos4[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
        if vec==R1 and float(s[2]) != float(0):
            pos1=[float(s[1])+float(shift1),float(s[2]),float(s[3])]
            for i in range(L):
                if i == 0:
                    r1=float(pos1[0])+shift*float(pos1[0])
                    r2=float(pos1[1])+shift*float(pos1[1])
                    r3=float(pos1[2])+shift*float(pos1[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos1[0])+shift*(i+1)*float(pos1[0])
                    r2=float(pos1[1])+shift*(i+1)*float(pos1[1])
                    r3=float(pos1[2])+shift*(i+1)*float(pos1[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
            pos2=[float(s[1])-float(shift1),float(s[2]),float(s[3])]
            for i in range(L):
                if i == 0:
                    r1=float(pos2[0])+shift*float(pos2[0])
                    r2=float(pos2[1])+shift*float(pos2[1])
                    r3=float(pos2[2])+shift*float(pos2[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos2[0])+shift*(i+1)*float(pos2[0])
                    r2=float(pos2[1])+shift*(i+1)*float(pos2[1])
                    r3=float(pos2[2])+shift*(i+1)*float(pos2[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
            pos3=[float(s[1]),float(s[2]),float(s[3])+float(shift1)]
            for i in range(L):
                if i == 0:
                    r1=float(pos3[0])+shift*float(pos3[0])
                    r2=float(pos3[1])+shift*float(pos3[1])
                    r3=float(pos3[2])+shift*float(pos3[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos3[0])+shift*(i+1)*float(pos3[0])
                    r2=float(pos3[1])+shift*(i+1)*float(pos3[1])
                    r3=float(pos3[2])+shift*(i+1)*float(pos3[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
            pos4=[float(s[1]),float(s[2]),float(s[3])-float(shift1)]
            for i in range(L):
                if i == 0:
                    r1=float(pos4[0])+shift*float(pos4[0])
                    r2=float(pos4[1])+shift*float(pos4[1])
                    r3=float(pos4[2])+shift*float(pos4[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos4[0])+shift*(i+1)*float(pos4[0])
                    r2=float(pos4[1])+shift*(i+1)*float(pos4[1])
                    r3=float(pos4[2])+shift*(i+1)*float(pos4[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
        if vec == R1 and float(s[3]) != float(0):
            pos1=[float(s[1])+float(shift1),float(s[2]),float(s[3])]
            for i in range(L):
                if i == 0:
                    r1=float(pos1[0])+shift*float(pos1[0])
                    r2=float(pos1[1])+shift*float(pos1[1])
                    r3=float(pos1[2])+shift*float(pos1[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos1[0])+shift*(i+1)*float(pos1[0])
                    r2=float(pos1[1])+shift*(i+1)*float(pos1[1])
                    r3=float(pos1[2])+shift*(i+1)*float(pos1[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
            pos2=[float(s[1])-float(shift1),float(s[2]),float(s[3])]
            for i in range(L):
                if i == 0:
                    r1=float(pos2[0])+shift*float(pos2[0])
                    r2=float(pos2[1])+shift*float(pos2[1])
                    r3=float(pos2[2])+shift*float(pos2[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos2[0])+shift*(i+1)*float(pos2[0])
                    r2=float(pos2[1])+shift*(i+1)*float(pos2[1])
                    r3=float(pos2[2])+shift*(i+1)*float(pos2[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
            pos3=[float(s[1]),float(s[2])+float(shift1),float(s[3])]
            for i in range(L):
                if i == 0:
                    r1=float(pos3[0])+shift*float(pos3[0])
                    r2=float(pos3[1])+shift*float(pos3[1])
                    r3=float(pos3[2])+shift*float(pos3[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos3[0])+shift*(i+1)*float(pos3[0])
                    r2=float(pos3[1])+shift*(i+1)*float(pos3[1])
                    r3=float(pos3[2])+shift*(i+1)*float(pos3[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))
            pos4=[float(s[1]),float(s[2])-float(shift1),float(s[3])]
            for i in range(L):
                if i == 0:
                    r1=float(pos4[0])+shift*float(pos4[0])
                    r2=float(pos4[1])+shift*float(pos4[1])
                    r3=float(pos4[2])+shift*float(pos4[2])
                    fout.write(("S %f %f %f\n")%(r1,r2,r3))
                else:
                    r1=float(pos4[0])+shift*(i+1)*float(pos4[0])
                    r2=float(pos4[1])+shift*(i+1)*float(pos4[1])
                    r3=float(pos4[2])+shift*(i+1)*float(pos4[2])
                    fout.write(("O %f %f %f\n")%(r1,r2,r3))

    
    
    
    
    
