from __future__ import division
import os
import sys
from math import sqrt
import numpy as np
import xml.etree.ElementTree as ET

def reorder_mats_xml(save,read):
	### INPUT FILE PARSE ###	 
	tree=ET.parse(read)
	root=tree.getroot()		
	### GET NODE INFORMATION ###	
	for child in root:
		### FIND ALL NODES ###
		ranks=child.findall('./')

		### GET V and Au INDICES FOR REORDERING ###
		rank=child.find('type').text
		rank=rank.splitlines()
		if(rank[1]=='V1'):
			print('#########################\nfirst V position is correct\n#########################')
			switch=False
		else:
			print('#########################\nneed to switch V positions\n#########################')
			switch=True

		rank=child.find('type').text
		rank=rank.splitlines()
		if(switch):
			print('switching V positions')
			counter=0
			indices=[]
      			indAu=[]
			firstAu=False
            		for i in range(len(rank)):
                        	if(rank[i].startswith('V')):
                                	counter+=1
                                	indices.append(i)
                        	if(rank[i]=='Au' and not(firstAu)):
                                	firstAu=True
                                	indAu.append(i)	
                        	if(rank[i]=='S'):
                                	firstAu=False

			print('\nthere are %s particles'% counter)
                	print('indices to switch are')
                	print(indices,indAu)
			
			tags=[]
			for node in ranks:
				tags.append(node.tag)
				
			for t in range(len(tags)):
				if(tags[t]!='box'):
					temp=child.find(tags[t]).text
					temp=temp.splitlines()
					if(len(temp)==len(rank)):
						temp[indices[0]],temp[indAu[0]]=temp[indAu[0]],temp[indices[0]]
						temp[indices[1]],temp[indAu[1]]=temp[indAu[1]],temp[indices[1]]
						#print(temp[indAu[0]],temp[indAu[1]])
						temp2='\n'.join(temp)
						newtemp=child.find(tags[t])
						newtemp.text=str(temp2)+'\n'
	########## WRITE THE NEW INDICES TO THE OUTPUT FILE ############
	tree.write(save)

def renumber_bodies_xml(save,read):
	### INPUT FILE PARSE ###	 
	tree=ET.parse(read)
	root=tree.getroot()		
	### GET NODE INFORMATION ###	
	for child in root:
		### GET V INDICES FOR RENAMING ###
		rank3=child.find('type').text
		rank3=rank3.splitlines()
		indices=[]
            	for i in range(len(rank3)):
                        if(rank3[i].startswith('V')):
                                indices.append(i)
				print(indices)

		rank=child.find('body').text
		rank=rank.splitlines()
		for i in range(len(rank)):
			if rank[i]=='1':
				rank[i]=str(indices[1]-1)
			else:
				pass
		rank2='\n'.join(rank)
		newrank=child.find('body')
		newrank.text=str(rank2)+'\n'

	########## WRITE THE NEW INDICES TO THE OUTPUT FILE ############
	tree.write(save)

###########################################################################
### This function simply copies an xml file, and renames the particle types
###########################################################################

def rename_particles_xml(save,read,newname,oldname):
    
    fid = open(save,'w')
    inputfile = open(read,'r')
    fdata=inputfile.read()
    data=fdata.splitlines()

    for line in data:
        s=line.split()
        if s[0]==oldname:    
            fid.write(('%s \n')%(newname))
        else:
            fid.write(('%s \n')%(line))

##################################################################
### This function simply copies an xml file, and renames the bonds
##################################################################

def rename_bonds_xml(save,read,newname,oldname):
    
    fid = open(save,'w')
    inputfile = open(read,'r')
    fdata=inputfile.read()
    data=fdata.splitlines()

    for line in data:
        if line.startswith(oldname):
            s=line.split()
            fid.write(('%s %s %s \n')%(newname,s[1],s[2]))
        else:
            fid.write(('%s \n')%(line))
            
            
######################################################################
### This function simply copies an xml file, and renames the dihedrals
######################################################################

def rename_dihedrals_xml(save,read,newname,oldname):
    
    fid = open(save,'w')
    inputfile = open(read,'r')
    fdata=inputfile.read()
    data=fdata.splitlines()

    for line in data:
        if line.startswith(oldname):
            s=line.split()
            fid.write(('%s %s %s %s\n')%(newname,s[1],s[2],s[3]))
        else:
            fid.write(('%s \n')%(line))

###################################################################
### This function simply copies an xml file, and renames the angles
###################################################################

def rename_angles_xml(save,read,newname,oldname):
    
    fid = open(save,'w')
    inputfile = open(read,'r')
    fdata=inputfile.read()
    data=fdata.splitlines()

    for line in data:
        if line.startswith(oldname):
            s=line.split()
            fid.write(('%s %s %s %s \n')%(newname,s[1],s[2],s[3]))
        else:
            fid.write(('%s \n')%(line))
            
##############################################################################
### This function shifts the positions of xyz files in the x y or z directions
##############################################################################

def shift_xml(save,read,shiftx,shifty,shiftz):


    #Open input file for reading in xml format

    fin=open(read,'r')
    fdata=fin.read()
    data=fdata.splitlines()

    #Open output file for writing in xml format

    fout=open(save,'w')
    
    writer=True
    for line in data:
        if writer==True:
            fout.write(('%s \n')%(line))
        if line.startswith('<position'):
            writer=False
        if writer==False and not(line.startswith('<position') or line.startswith('</position')):              
            s=line.split()
            x=float(s[0])+shiftx
            y=float(s[1])+shifty
            z=float(s[2])+shiftz
            fout.write(("%f %f %f\n")%(x,y,z))
        if line.startswith('</position'):
            writer=True
            
##############################################################################
### This function scales the positions of xyz files in the x y or z directions
##############################################################################           

def scale_xml(save,read,scalex,scaley,scalez):


    #Open input file for reading in xml format

    fin=open(read,'r')
    fdata=fin.read()
    data=fdata.splitlines()

    #Open output file for writing in xml format

    fout=open(save,'w')
    
    writer=True
    for line in data:
        if writer==True:
            fout.write(('%s \n')%(line))
        if line.startswith('<position'):
            writer=False
        if writer==False and not(line.startswith('<position') or line.startswith('</position')):              
            s=line.split()
            x=float(s[0])*scalex
            y=float(s[1])*scaley
            z=float(s[2])*scalez
            fout.write(("%f %f %f\n")%(x,y,z))
        if line.startswith('</position'):
            writer=True
            
##########################
### Example run script ###
##########################
"""        
import os
import sys

import manipulate_xml as m

save='newname.xml'
read='tetra.xml'
newname='Au'
oldname='C'

m.rename_particles_xml(save,read,newname,oldname)
"""
