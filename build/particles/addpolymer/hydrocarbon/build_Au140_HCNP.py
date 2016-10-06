from __future__ import division
import numpy as np
from numpy import linalg as LA
import rotate_polymer as rp

def build_HC_chain_noS(save,n):
    
    fid = open(save,'w')
    
    CH2=[0.0,0.0,0.0]
    S=[-1.58237,0.899165,0.0]
    Spos=np.add(CH2,S)
    CH2vec1=[1.24946,0.883031,0.0]
    CH2vec2=[1.24946,-0.883031,0.0]
    
    fid.write(str(1+n)+'\n\n')
    current= Spos
    for  x in range(0,n):
        if(x==0):
            current=CH2
        if(x%2==0):
            current=np.add(current,CH2vec1)
        else:
            current=np.add(current,CH2vec2)
        fid.write(('CH2 %s %s %s\n')%(current[0],current[1],current[2]))    
    if(n%2 == 0):
        newcurrent=np.add(current,CH2vec1)
        fid.write(('CH3 %s %s %s\n')%(newcurrent[0],newcurrent[1],newcurrent[2]))
    else:
        newcurrent=np.add(current,CH2vec2)
        fid.write(('CH3 %s %s %s\n')%(newcurrent[0],newcurrent[1],newcurrent[2]))


def graft_HC_hex(save,readpoly,readnp,Rcf,Rct):
    
    fin1 = open(readpoly,'r')
    fin2 = open(readnp,'r')
    fout = open(save,'w')
    
    
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[2:]
    print('length of poly is '+str(len(data1)+1))
    
    f2data=fin2.read()
    data2=f2data.splitlines()
    data2=data2[2:]

    centers=np.array([[0.0,0.0,0.0]])
    centersout=np.array([[0.0,0.0,0.0]])
    tophexes=np.array([[0.0,0.0,0.0]])
    bothexes=np.array([[0.0,0.0,0.0]])
    
    count=int(8*len(data1)+1)        
    centers=np.append(centers,[[Rcf,Rcf,Rcf]],axis=0)    
    centers=np.append(centers,[[Rcf,Rcf,-Rcf]],axis=0)
    centers=np.append(centers,[[Rcf,-Rcf,-Rcf]],axis=0)
    centers=np.append(centers,[[Rcf,-Rcf,Rcf]],axis=0)    
    centers=np.append(centers,[[-Rcf,Rcf,Rcf]],axis=0)
    centers=np.append(centers,[[-Rcf,-Rcf,Rcf]],axis=0)
    centers=np.append(centers,[[-Rcf,Rcf,-Rcf]],axis=0)
    centers=np.append(centers,[[-Rcf,-Rcf,-Rcf]],axis=0)
    
    out=1.137380030    
    centersout=np.append(centersout,[[Rcf+out,Rcf+out,Rcf+out]],axis=0)    
    centersout=np.append(centersout,[[Rcf+out,Rcf+out,-Rcf-out]],axis=0)
    centersout=np.append(centersout,[[Rcf+out,-Rcf-out,-Rcf-out]],axis=0)
    centersout=np.append(centersout,[[Rcf+out,-Rcf-out,Rcf+out]],axis=0)    
    centersout=np.append(centersout,[[-Rcf-out,Rcf+out,Rcf+out]],axis=0)
    centersout=np.append(centersout,[[-Rcf-out,-Rcf-out,Rcf+out]],axis=0)
    centersout=np.append(centersout,[[-Rcf-out,Rcf+out,-Rcf-out]],axis=0)
    centersout=np.append(centersout,[[-Rcf-out,-Rcf-out,-Rcf-out]],axis=0)
    
    
    for line in data2:
        s=line.split()
        csv=[s[1],s[2],s[3]]
        cf=LA.norm(csv)
        vec=np.sqrt((float(s[1])**2)+(float(s[2])**2)+(float(s[3])**2))
        if vec == Rct and float(s[3]) >= 0.0 and float(s[2]) >= 0 and float(s[1]) >= 0:
            count+=int((2/3.)*len(data1)+1)
            r1=float(s[1])
            r2=float(s[2])
            r3=float(s[3])
            v=[[r1,r2,r3]]
            for rv in range(1,centers.shape[0]):
                diff=np.subtract(centers[rv],v)
                csv=[centers[rv][0],centers[rv][1],centers[rv][2]]
                cf=LA.norm(csv)
                if LA.norm(diff)<=2.0:
                    cfv=[float(csv[0])/float(cf),float(csv[1])/float(cf),float(csv[2])/float(cf)]
                    addout=np.multiply(1.97,cfv)
                    tophexes=np.append(tophexes,np.add(addout,np.subtract(v,1.393*diff)),axis=0)
                    tophexes=np.append(tophexes,np.add(addout,np.add(v,4*diff)),axis=0)
        if vec == Rct and float(s[3]) >= 0.0 and float(s[2]) <= 0 and float(s[1]) <= 0:
            count+=int((2/3.)*len(data1)+1)
            r1=float(s[1])
            r2=float(s[2])
            r3=float(s[3])
            v=[[r1,r2,r3]]
            for rv in range(1,centers.shape[0]):
                diff=np.subtract(centers[rv],v)                
                csv=[centers[rv][0],centers[rv][1],centers[rv][2]]
                cf=LA.norm(csv)
                if LA.norm(diff)<=2.0:
                    cfv=[float(csv[0])/float(cf),float(csv[1])/float(cf),float(csv[2])/float(cf)]
                    addout=np.multiply(1.97,cfv)
                    tophexes=np.append(tophexes,np.add(addout,np.subtract(v,1.393*diff)),axis=0)
                    tophexes=np.append(tophexes,np.add(addout,np.add(v,4*diff)),axis=0)
        if vec == Rct and float(s[3]) <= 0.0 and float(s[2]) >= 0 and float(s[1]) <= 0:
            count+=int((2/3.)*len(data1)+1)
            r1=float(s[1])
            r2=float(s[2])
            r3=float(s[3])
            v=[[r1,r2,r3]]
            for rv in range(1,centers.shape[0]):
                diff=np.subtract(centers[rv],v)                
                csv=[centers[rv][0],centers[rv][1],centers[rv][2]]
                cf=LA.norm(csv)
                if LA.norm(diff)<=2.0:
                    cfv=[float(csv[0])/float(cf),float(csv[1])/float(cf),float(csv[2])/float(cf)]
                    addout=np.multiply(1.97,cfv)
                    tophexes=np.append(tophexes,np.add(addout,np.subtract(v,1.393*diff)),axis=0)
                    tophexes=np.append(tophexes,np.add(addout,np.add(v,4*diff)),axis=0)
        if vec == Rct and float(s[3]) <= 0.0 and float(s[2]) <= 0 and float(s[1]) >= 0:
            count+=int((2/3.)*len(data1)+1)
            r1=float(s[1])
            r2=float(s[2])
            r3=float(s[3])
            v=[[r1,r2,r3]]
            for rv in range(1,centers.shape[0]):
                diff=np.subtract(centers[rv],v)                
                csv=[centers[rv][0],centers[rv][1],centers[rv][2]]
                cf=LA.norm(csv)
                if LA.norm(diff)<=2.0:
                    cfv=[float(csv[0])/float(cf),float(csv[1])/float(cf),float(csv[2])/float(cf)]
                    addout=np.multiply(1.97,cfv)
                    tophexes=np.append(tophexes,np.add(addout,np.subtract(v,1.393*diff)),axis=0)
                    tophexes=np.append(tophexes,np.add(addout,np.add(v,4*diff)),axis=0)
        V=[[0,0,0]]
        if vec == Rct and float(s[3]) >= 0.0 and float(s[2]) >= 0 and float(s[1]) <= 0:
            count+=int((1/3.)*len(data1)+1)
            r1=float(s[1])
            r2=float(s[2])
            r3=float(s[3])
            v=[[r1,r2,r3]]
            for rv in range(1,centers.shape[0]):
                k=np.subtract(V,centers[rv])
                knor=LA.norm(k)
                kt=[float(k[0][0])/float(knor),float(k[0][1])/float(knor),float(k[0][2])/float(knor)]               
                K=[[0,-kt[2],kt[1]],[kt[2],0,-kt[0]],[-kt[1],kt[0],0]]
                ###                
                csv=[centers[rv][0],centers[rv][1],centers[rv][2]]
                cf=LA.norm(csv)
                cfv=[float(csv[0])/float(cf),float(csv[1])/float(cf),float(csv[2])/float(cf)]
                ###
                v1=np.subtract(centers[rv],v)
                vn=[v1[0][0],v1[0][1],v1[0][2]]
                vnor=LA.norm(vn)
                vt=[float(vn[0])/float(vnor),float(vn[1])/float(vnor),float(vn[2])/float(vnor)]
                ###
                theta=np.pi/7
                theta2=2*np.pi-theta
                R=np.add(np.identity(3),np.add(np.multiply(np.sin(theta),K),np.multiply((1-np.cos(theta)),np.dot(K,K))))
                R2=np.add(np.identity(3),np.add(np.multiply(np.sin(theta2),K),np.multiply((1-np.cos(theta2)),np.dot(K,K))))
                if LA.norm(v1)<=2.0:
                    addout=np.multiply(1.97,cfv)                     
                    nnfac=4.49429  
                    u1=np.add(np.multiply(Rct,cfv),np.multiply(nnfac,vt))
                    u=[u1[0],u1[1],u1[2]]
                    vr=np.dot(R,u)
                    vr2=np.dot(R2,u)
                    bothexes=np.append(bothexes,np.add([addout],[vr]),axis=0)
                    bothexes=np.append(bothexes,np.add([addout],[vr2]),axis=0)  
        V=[[0,0,0]]
        if vec == Rct and float(s[3]) >= 0.0 and float(s[2]) <= 0 and float(s[1]) >= 0:
            count+=int((1/3.)*len(data1)+1)
            r1=float(s[1])
            r2=float(s[2])
            r3=float(s[3])
            v=[[r1,r2,r3]]
            for rv in range(1,centers.shape[0]):
                k=np.subtract(V,centers[rv])
                knor=LA.norm(k)
                kt=[float(k[0][0])/float(knor),float(k[0][1])/float(knor),float(k[0][2])/float(knor)]               
                K=[[0,-kt[2],kt[1]],[kt[2],0,-kt[0]],[-kt[1],kt[0],0]]
                ###                
                csv=[centers[rv][0],centers[rv][1],centers[rv][2]]
                cf=LA.norm(csv)
                cfv=[float(csv[0])/float(cf),float(csv[1])/float(cf),float(csv[2])/float(cf)]
                ###
                v1=np.subtract(centers[rv],v)
                vn=[v1[0][0],v1[0][1],v1[0][2]]
                vnor=LA.norm(vn)
                vt=[float(vn[0])/float(vnor),float(vn[1])/float(vnor),float(vn[2])/float(vnor)]
                ###
                theta=np.pi/7
                theta2=2*np.pi-theta
                R=np.add(np.identity(3),np.add(np.multiply(np.sin(theta),K),np.multiply((1-np.cos(theta)),np.dot(K,K))))
                R2=np.add(np.identity(3),np.add(np.multiply(np.sin(theta2),K),np.multiply((1-np.cos(theta2)),np.dot(K,K))))
                if LA.norm(v1)<=2.0:
                    addout=np.multiply(1.97,cfv)                     
                    nnfac=4.49429  
                    u1=np.add(np.multiply(Rct,cfv),np.multiply(nnfac,vt))
                    u=[u1[0],u1[1],u1[2]]
                    vr=np.dot(R,u)
                    vr2=np.dot(R2,u)
                    bothexes=np.append(bothexes,np.add([addout],[vr]),axis=0)
                    bothexes=np.append(bothexes,np.add([addout],[vr2]),axis=0)  
        V=[[0,0,0]]
        if vec == Rct and float(s[3]) <= 0.0 and float(s[2]) <= 0 and float(s[1]) <= 0:
            count+=int((1/3.)*len(data1)+1)
            r1=float(s[1])
            r2=float(s[2])
            r3=float(s[3])
            v=[[r1,r2,r3]]
            for rv in range(1,centers.shape[0]):
                k=np.subtract(V,centers[rv])
                knor=LA.norm(k)
                kt=[float(k[0][0])/float(knor),float(k[0][1])/float(knor),float(k[0][2])/float(knor)]               
                K=[[0,-kt[2],kt[1]],[kt[2],0,-kt[0]],[-kt[1],kt[0],0]]
                ###                
                csv=[centers[rv][0],centers[rv][1],centers[rv][2]]
                cf=LA.norm(csv)
                cfv=[float(csv[0])/float(cf),float(csv[1])/float(cf),float(csv[2])/float(cf)]
                ###
                v1=np.subtract(centers[rv],v)
                vn=[v1[0][0],v1[0][1],v1[0][2]]
                vnor=LA.norm(vn)
                vt=[float(vn[0])/float(vnor),float(vn[1])/float(vnor),float(vn[2])/float(vnor)]
                ###
                theta=np.pi/7
                theta2=2*np.pi-theta
                R=np.add(np.identity(3),np.add(np.multiply(np.sin(theta),K),np.multiply((1-np.cos(theta)),np.dot(K,K))))
                R2=np.add(np.identity(3),np.add(np.multiply(np.sin(theta2),K),np.multiply((1-np.cos(theta2)),np.dot(K,K))))
                if LA.norm(v1)<=2.0:
                    addout=np.multiply(1.97,cfv)                     
                    nnfac=4.49429  
                    u1=np.add(np.multiply(Rct,cfv),np.multiply(nnfac,vt))
                    u=[u1[0],u1[1],u1[2]]
                    vr=np.dot(R,u)
                    vr2=np.dot(R2,u)
                    bothexes=np.append(bothexes,np.add([addout],[vr]),axis=0)
                    bothexes=np.append(bothexes,np.add([addout],[vr2]),axis=0)    
        V=[[0,0,0]]
        if vec == Rct and float(s[3]) <= 0.0 and float(s[2]) >= 0 and float(s[1]) >= 0:
            count+=int((1/3.)*len(data1)+1)
            r1=float(s[1])
            r2=float(s[2])
            r3=float(s[3])
            v=[[r1,r2,r3]]
            for rv in range(1,centers.shape[0]):
                k=np.subtract(V,centers[rv])
                knor=LA.norm(k)
                kt=[float(k[0][0])/float(knor),float(k[0][1])/float(knor),float(k[0][2])/float(knor)]               
                K=[[0,-kt[2],kt[1]],[kt[2],0,-kt[0]],[-kt[1],kt[0],0]]
                ###                
                csv=[centers[rv][0],centers[rv][1],centers[rv][2]]
                cf=LA.norm(csv)
                cfv=[float(csv[0])/float(cf),float(csv[1])/float(cf),float(csv[2])/float(cf)]
                ###
                v1=np.subtract(centers[rv],v)
                vn=[v1[0][0],v1[0][1],v1[0][2]]
                vnor=LA.norm(vn)
                vt=[float(vn[0])/float(vnor),float(vn[1])/float(vnor),float(vn[2])/float(vnor)]
                ###
                theta=np.pi/7
                theta2=2*np.pi-theta
                R=np.add(np.identity(3),np.add(np.multiply(np.sin(theta),K),np.multiply((1-np.cos(theta)),np.dot(K,K))))
                R2=np.add(np.identity(3),np.add(np.multiply(np.sin(theta2),K),np.multiply((1-np.cos(theta2)),np.dot(K,K))))
                if LA.norm(v1)<=2.0:
                    addout=np.multiply(1.97,cfv)                     
                    nnfac=4.49429  
                    u1=np.add(np.multiply(Rct,cfv),np.multiply(nnfac,vt))
                    u=[u1[0],u1[1],u1[2]]
                    vr=np.dot(R,u)
                    vr2=np.dot(R2,u)
                    bothexes=np.append(bothexes,np.add([addout],[vr]),axis=0)
                    bothexes=np.append(bothexes,np.add([addout],[vr2]),axis=0)    
                    
    fout.write(('%s\n\n')%(str(len(data2)+count)))
    
    for line in data2:
        s=line.split()
        fout.write(('%s %f %f %f\n')%(s[0],float(s[1]),float(s[2]),float(s[3])))
    
    for v in range(1,centersout.shape[0]):
        P=rp.align_vector(readpoly,centersout[v])
        for i in range(P.shape[0]):
            P[i]=np.add(centersout[v],P[i])
            if(i==P.shape[0]-1):
                fout.write(('CH3 %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
            elif(i==0):
                fout.write(('S %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
            else:
                fout.write(('CH2 %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
                
    for v in range(1,tophexes.shape[0]):
        P=rp.align_vector(readpoly,tophexes[v])
        for i in range(P.shape[0]):
            P[i]=np.add(tophexes[v],P[i])
            if(i==P.shape[0]-1):
                fout.write(('CH3 %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
            elif(i==0):
                fout.write(('S %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
            else:
                fout.write(('CH2 %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))              

    for v in range(1,bothexes.shape[0]):
        P=rp.align_vector(readpoly,bothexes[v])
        for i in range(P.shape[0]):
            P[i]=np.add(bothexes[v],P[i])
            if(i==P.shape[0]-1):
                fout.write(('CH3 %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
            elif(i==0):
                fout.write(('S %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))
            else:
                fout.write(('CH2 %f %f %f\n')%(P[i][0],P[i][1],P[i][2]))    

def graft_square_faces(save,readpoly,readnp,R_face):
    fin1 = open(readpoly,'r')
    fin2 = open(readnp,'r')
    fout = open(save,'w')
    
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[2:]
    
    f2data=fin2.read()
    data2=f2data.splitlines()
    data2=data2[2:]
    #####################
    #AuS is length radial outward from face for Au-S distance
    #####################
    AuS=1.97
    
    vecs=np.array([[0.0,0.0,0.0]])
    count=(6*len(data1)+6)
    
    Spos=R_face+AuS
    vecs=np.append(vecs,[[0.0,0.0,Spos]],axis=0)
    vecs=np.append(vecs,[[0.0,0.0,-Spos]],axis=0)
    vecs=np.append(vecs,[[0.0,Spos,0.0]],axis=0)
    vecs=np.append(vecs,[[0.0,-Spos,0.0]],axis=0) 
    vecs=np.append(vecs,[[Spos,0.0,0.0]],axis=0)
    vecs=np.append(vecs,[[-Spos,0.0,0.0]],axis=0)
    
    fout.write(('%s\n\n')%(str(len(data2)+count)))
    
    for line in data2:
        s=line.split()
        fout.write(('%s %f %f %f\n')%(s[0],float(s[1]),float(s[2]),float(s[3])))
    
    for v in range(1,vecs.shape[0]):
        P=rp.align_vector(readpoly,vecs[v])
        for x in range(P.shape[0]):
            P[x]=np.add(vecs[v],P[x])
            if(x==P.shape[0]-1):
                fout.write(('CH3 %f %f %f\n')%(P[x][0],P[x][1],P[x][2]))
            elif(x==0):
                fout.write(('S %f %f %f\n')%(P[x][0],P[x][1],P[x][2]))
            else:
                fout.write(('CH2 %f %f %f\n')%(P[x][0],P[x][1],P[x][2]))


