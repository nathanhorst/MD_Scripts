# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 16:08:55 2016

@author: waltmann
"""


import numpy as np

def get_energies(mylog):
    with open(mylog) as f:
        data = f.readlines()
        data=data[2:]
    energy=np.array([0.0])    
    for x in range(0,len(data)):
        energy=np.append(energy,[float(data[x].split()[2])*0.002872],axis=0)
    return energy[1:]


def get_positions(Hfile):
	with open(Hfile) as f:
		data=f.readlines()
	pos = np.array([0.0])
	for i in range(len(data)):
		pos=np.append(pos,[float(data[i])],axis=0)
	return pos[1:]

def correlation(numks, energy):
    correlation=np.array([[0.0,1.0]])    
    for k in range(0,numks):
        #print k
        t=energy.shape[0]-k
        etk=0
        et2total=0
        ettotal=0
	corr=0
        for i in range(0,t):
            et2total+=(energy[i]**2)
            ettotal+=energy[i]
            etk+=energy[i]*energy[i+k]
        corr=((etk/t)-(ettotal/t)**2)/(et2total/t-(ettotal/t)**2)
        correlation=np.append(correlation,[[k,corr]],axis=0)
    return correlation

def estimated_autocorrelation(x):
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(np.arange(n, 0, -1)))
    return result
    
    
def revised_correlation(numks,energy):
    correlation=np.array([[0.0,1.0]]) 
    for k in range(0,numks):
        print k
        top=0
        avg=np.average(energy[:k])
        if(k==0):
            avg=np.average(energy)
        #print avg
        t=energy.shape[0]-k
	print energy.shape[0]
        for i in range(0,t):
            top+=(energy[i]*energy[i+k]-avg**2)/(energy[i]**2-avg**2)
        correlation=np.append(correlation,[[k,top/t]],axis=0)
    return correlation[1:]
    
