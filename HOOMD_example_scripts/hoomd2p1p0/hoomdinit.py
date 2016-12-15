from __future__ import division
import numpy as np


def momentI(rigidposmatrix,rigidmassmatrix):
	P=rigidposmatrix
	M=rigidmassmatrix
	### GET CENTER OF MASS ###
	Mtot=sum(M)
	print 'Mtot\n',Mtot
	cx,cy,cz=[0.0],[0.0],[0.0]
	for i in range(len(M)):
		cx=np.append(cx,M[i]*P[i][0])
		cy=np.append(cy,M[i]*P[i][1])		
		cz=np.append(cz,M[i]*P[i][2])
	cX = (sum(cx))/Mtot
	cY = (sum(cy))/Mtot
	cZ = (sum(cz))/Mtot
	com=[cX,cY,cZ]
	print 'com\n',com
	### GET MOMENTS OF INERTIA ###
	cXYZ=[[0.0,0.0,0.0]]
	for i in range(len(M)):
		cXYZ = np.append(cXYZ,[[P[i][0]-cX,P[i][1]-cY,P[i][2]-cZ]],axis=0)
		#cXYZ = np.append(cXYZ,[[P[i][0],P[i][1],P[i][2]]],axis=0)
	cXYZ=cXYZ[1:]
	print 'cXYZ\n',cXYZ
	Ixx,Iyy,Izz = [0.0],[0.0],[0.0]
	for i in xrange(len(M)):
		Ixx = np.append(Ixx,M[i]*(cXYZ[i][1]**2 + cXYZ[i][2]**2))
		Iyy = np.append(Iyy,M[i]*(cXYZ[i][0]**2 + cXYZ[i][2]**2))
		Izz = np.append(Izz,M[i]*(cXYZ[i][0]**2 + cXYZ[i][1]**2))
		Imatrix = np.array([sum(Ixx),sum(Iyy),sum(Izz)])

	print 'Imatrix\n',Imatrix,'\nCoM\n',com
	return Imatrix,com

"""
returns  a matrix containing the posiitons of all the V particle in the order
they are found in the file.
"""

def v_pos_matrix(inputfile):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[4:]     
    
    tipdone=False
    posdone=False
    pos=False
    tip=False
    tipcount=0
    poscount=0
    vnums=np.array([-5])
    position= np.array([[0.0,0.0,0.0]])
    v_position= np.array([[0.0,0.0,0.0]])
    for line in data1:
        #print line
        if(tip):
            tipcount+=1
        if(pos):
            poscount+=1
        s=line.split()
        if(s[0]=='</type>'):
            tip=False
            tipdone=True
        elif(tip==True and s[0][:1]=='V'):
            vnums=np.append(vnums,[tipcount],axis=0)
        elif(s[0][:5]=='<type'):
            tip=True
        elif(s[0]=='</position>'):
            pos=False 
            posdone=True
        elif(pos):
            vec=[float(s[0]),float(s[1]),float(s[2])]
            position=np.append(position,[vec], axis=0)
        elif(s[0][:9]=='<position'):
            pos=True
        elif(posdone and tipdone):
            for x in range(1,len(vnums)):
                v_position=np.append(v_position,[position[vnums[x]]], axis=0)
            #print v_position    
            return v_position

###10000 particle max
def type_pos_matrix_by_particle(inputfile,typ):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[4:]     
    v=v_pos_matrix(inputfile)
    v=v[1:]
    tipdone=False
    posdone=False
    pos=False
    tip=False
    tipcount=0
    poscount=0
    nums=np.array([-5])
    vnums=np.array([0])
    vcount=0
    position= np.array([[0.0,0.0,0.0]])
    v_position= np.array([[0.0,0.0,0.0]])
    for line in data1:
        #print line
        if(tip):
            tipcount+=1
        if(pos):
            poscount+=1
        s=line.split()
        if(s[0]=='</type>'):
	    vnums=np.append(vnums,[vcount],axis=0)
	    print vcount  
            vcount=0 
            tip=False
            tipdone=True
        elif(tip==True and s[0]==typ):
            nums=np.append(nums,[tipcount],axis=0)
            vcount+=1
        elif(tip==True and s[0][:1]=='V'):
            vnums=np.append(vnums,[vcount],axis=0)
	    print vcount  
            vcount=0
        elif(s[0][:5]=='<type'):
            tip=True
        elif(s[0]=='</position>' ):
            pos=False 
            posdone=True
        elif(pos):
            vec=[float(s[0]),float(s[1]),float(s[2])]
            position=np.append(position,[vec], axis=0)
        elif(s[0][:9]=='<position'):
            pos=True
        elif(posdone and tipdone):
            for x in range(1,len(nums)):
                v_position=np.append(v_position,[position[nums[x]]], axis=0)
    	    vnums=vnums[2:]
    	    v_position=v_position[1:]
    	    vs=np.zeros((len(v),np.amax(vnums),3))
            vposcount=0
    	    for i in range(0,len(vs)):
        	for x in range(0,vnums[i]):
            		vs[i][x]=v_position[vposcount]
            		vposcount+=1
    	    return vs,vnums


def set_up_rigid(system,sysfile,rigid,distance_scale=1, massAu=14.0424, massS=2.28601):
    tags = []
    print distance_scale
    au,aunums=type_pos_matrix_by_particle(sysfile,'Au')
    s,snums=type_pos_matrix_by_particle(sysfile,'S')
    v=v_pos_matrix(sysfile)
    v=v[1:]
    for b in system.particles:
        if b.type[:1] == 'V':
            tags.append(b.tag)

    for i in range(0,len(aunums)):
        x=np.array([[0,0,0]])
        for tr in range(0,aunums[i]):
            x=np.append(x,[au[i][tr]],axis=0)
        x=x[1:]
        q=np.array([[0,0,0]])
        for tr in range(0,snums[i]):
            q=np.append(q,[s[i][tr]],axis=0)
        q=q[1:]
        g=np.append(x,q,axis=0)
        t=['Au']
        for poo in range(len(x)-1): 
            t.append('Au')
        for b in range(0,len(q)):
    		t.append('S')
        g=np.subtract(g,v[i])
        aumass=np.full(aunums[i],massAu)
        smass=np.full(snums[i],massS)
        masses=np.concatenate((aumass,smass),axis=0)
        print("positions and masses")
        print g
        print  masses
        moment,com=momentI(g*distance_scale,masses)
        system.particles.get(tags[i]).moment_inertia=moment
        system.particles.get(tags[i]).velocity=[0.0,0.0,0.0]
        #system.particles.get(tags[i]).orientation=[1.0,0.0,0.0,0.0]
        rigid.set_param(str('V'+str(i+1)), positions=g, types=t)

def vector_length(x,y,z):
    c=np.sqrt(x**2+y**2+z**2)
    return c

"""
tkaes a vetor and returns its length
"""   
def length(vector):
    return vector_length(vector[0],vector[1],vector[2])   

"""
does the distance between particles without knowing the box size expliciitly.
just the file it came from
"""
def part_distance(part1,part2,inputfile):
    fin1= open(inputfile,'r')
    f1data=fin1.read()
    data1=f1data.splitlines()
    #print data1[3]
    s=data1[3].split()
    #print s
    lx=0
    ly=0
    lz=0
    xy=0
    for i in range(len(s)):
        if s[i][:2]=='lx':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            lx=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='ly':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            ly=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='lz':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            lz=float(s[i][count+1:b+count-len(s[i])])
        elif s[i][:2]=='xy':
            count=0
            while(s[i][count]!='"'):
                count+=1
            b=1
            while(s[i][count+b]!='"'):
                b+=1
            xy=float(s[i][count+1:b+count-len(s[i])])
    #print lx
    #print ly
    #print lz
    #print part1
    #print part2[0]
    x=abs(float(part1[0])-float(part2[0]))
    if(x>(lx/2)-1):
        x=lx-x
    y=abs(float(part1[1])-float(part2[1]))
    
    z=abs(float(part1[2])-float(part2[2]))
    if(z>(lz/2)-1):
        z=lz-z
    #print x
    #print y 
    #print z
    vector=[x,y,z]
    #print vector
    v2=[xy*ly,ly,0]
    v2hat=v2/length(v2)
    ycheck=np.dot((vector),v2hat)
    #print('v2hat is '+ str(v2hat))
    #print np.subtract(part1,part2)
    #print('ycheck is '+ str(ycheck))
    if(ycheck>(length(v2)/2)-1):
        vector=np.subtract(vector,v2)
        print vector
    if(vector[0]>(lx/2)-1):
        vector[0]=lx-vector[0]
        print vector
    return vector_length(vector[0],vector[1],vector[2])

 
