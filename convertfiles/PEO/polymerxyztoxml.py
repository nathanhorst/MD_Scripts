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
        
    def bonds(fid, data, Lpoly):
        fid.write('<bond>\n')
        i=-1
        for line in data:
            i += 1
            s = line.split()
            if s[0]=='S':
                 for j in range(1,Lpoly):                
                    fid.write(("fene %d %d\n")%(int(i+j-1),int(i+j)))
        fid.write('</bond>\n')   
                    
                
    particle_type(fid,data)
    position(fid,data)
    body(fid,data)
#    mass(fid,data)
    bonds(fid,data,Lpoly)  
    fid.write('</configuration>\n</hoomd_xml>')
    fid.close()
    inputfile.close()  
