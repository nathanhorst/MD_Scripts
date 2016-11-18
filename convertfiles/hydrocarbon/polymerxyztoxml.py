import os
import sys

import numpy as np

def polynp_conv(L,save,read,Lpoly):
    
    fid = open(save,'w')
    inputfile = open(read,'r')
    fdata=inputfile.read()
    
    data=fdata.splitlines()
    data=data[2:]
    
    # first lines in xml file
    fid.write('<?xml version="1.0" encoding="UTF-8"?>\n<hoomd_xml version="1.1">\n<configuration time_step="0">\n')
    #box size
    fid.write(('<box units="sigma"  lx="%f" ly="%f" lz="%f"/>\n')%(L[0],L[1],L[2]))

    def particle_type(fid, data):
        fid.write('<type>\n')
        for line in data:
            s = line.split()
            fid.write(('%s\n')%(s[0]))
        fid.write('</type>\n')   
        
    def position(fid, data):
        fid.write('<position units="sigma" >\n')
        for line in data:
            s = line.split()
            fid.write(("%s %s %s\n")%(s[1],s[2],s[3]))
        fid.write('</position>\n')
        
    def body(fid, data):
        fid.write('<body>\n')
        count=0
        canCount=False
        for line in data:
            s = line.split()
            if s[0] == 'Au' or s[0] == 'S' or s[0]=='V':
	    	if(s[0]=='Au' and canCount):
			count+=1
		canCount=False
		fid.write(str(count)+'\n')
            else:
                fid.write('-1\n')
		canCount=True
        fid.write('</body>\n')
        
    def mass(fid, data):
        fid.write('<mass>\n')
        for line in data:
            s = line.split()
            if s[0] == 'Au'or s[0]=='V':
                fid.write('14.0424\n')
            if s[0] == 'S':
                fid.write('2.28601\n')            
            if s[0] == 'CH2':
                fid.write('1\n')
            if s[0] == 'CH3':
                fid.write('1.07199\n')
        fid.write('</mass>\n')
        
    def bonds(fid, data, Lpoly,PFconfig=False):
        fid.write('<bond>\n')
        if PFconfig:
                lol=-1
                vpos=np.array([0.0])
                for line in data:
                        lol+=1
                        s=line.split()
                        if s[0] == 'V':
                                vpos=np.append(vpos,lol)
                fid.write('V-V '+str(int(vpos[1]))+' '+str(int(vpos[2]))+'\n')
        i=-1
        for line in data:
            i += 1
            s = line.split()
            if s[0]=='S':
                 for j in range(1,2):                
                    fid.write(("S-C %d %d\n")%(int(i+j-1),int(i+j)))
                 for j in range(2,Lpoly):                
                    fid.write(("C-C %d %d\n")%(int(i+j-1),int(i+j)))
        fid.write('</bond>\n')   

    def angles(fid, data, Lpoly):
        fid.write('<angle>\n')
        i=-1
        for line in data:
            i+=1
            s = line.split()
	    if s[0]=='V':
		vcoord=i
            if s[0]=='S':
                for j in range(0,Lpoly-2):
                    if j==0:
			#fid.write(("Au-S-CH2 %d %d %d\n")%(int(vcoord),int(i),int(i+1)))
                        fid.write(("S-CH2-CH2 %d %d %d\n")%(int(i),int(i+1),int(i+2)))
		    elif j==Lpoly-3:
			fid.write(("CH2-CH2-CH3 %d %d %d\n")%(int(i+j),int(i+j+1),int(i+j+2)))
                    else:
                        fid.write(("CH2-CH2-CH2 %d %d %d\n")%(int(i+j),int(i+j+1),int(i+j+2)))
        fid.write('</angle>\n')
        
    def dihedrals(fid, data, Lpoly):
        fid.write('<dihedral>\n')
        i=-1
        for line in data:
            i+=1
            s = line.split()
            if s[0]=='S':
                for j in range(0,Lpoly-3):
                    fid.write(("phi1 %d %d %d %d\n")%(int(i+j),int(i+j+1),int(i+j+2),int(i+j+3)))
        fid.write('</dihedral>\n')

    particle_type(fid,data)
    position(fid,data)
    body(fid,data)
    mass(fid,data)
    bonds(fid,data,Lpoly)
    angles(fid,data,Lpoly)
    dihedrals(fid,data,Lpoly)    
    fid.write('</configuration>\n</hoomd_xml>')
    fid.close()
    inputfile.close()  
