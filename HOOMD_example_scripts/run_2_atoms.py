"""
Created on Tue Jul  5 10:44:23 2016
@author: waltmann
"""

from __future__ import division
from hoomd_script import *
import math
import sys
sys.path.insert(0, "../analysis/hydrocarbon/")
print sys.path
print '\n\n\n\n'
import util as u
import V_distance as vd
import os
import numpy as np
kb=8.61733e-5*30
distance_scale=1.0/3.96
real_temp=30
fs=np.array([[0.0,0.0]])
timesteps=0

"""
dirac potential
"""
def dirac_potential(r,rmin,rmax,r0):
	a=.5
	V=300*(-1/(a * np.sqrt(np.pi)))* np.exp(-1*(r-r0)**2/a**2)
	F=300*(-2*(r-r0)/((a**3)* np.sqrt(np.pi)))*np.exp(-1*(r-r0)**2/a**2)
	return (V,F)
def lj(r, rmin, rmax, epsilon, sigma):
    V = 4 * epsilon * ( (sigma / r)**12 - (sigma / r)**6);
    F = 4 * epsilon / r * ( 12 * (sigma / r)**12 - 6 * (sigma / r)**6);
    return (V, F)
def dirac_bond(r,rmin,rmax,kappa,r0):
	a=.5
	V=kappa*(-1/(a * np.sqrt(np.pi)))* np.exp(-1*(r-r0)**2/a**2)
	F=kappa*(-2*(r-r0)/((a**3)* np.sqrt(np.pi)))*np.exp(-1*(r-r0)**2/a**2)
	return (V,F)
def harmonic(r, rmin, rmax, kappa, r0):
   V = 0.5 * kappa * (r-r0)**2;
   F = -kappa*(r-r0);
   return (V, F)
    
sysfile= "2nano.xml"
#sysfile='Au201shellHC.xml'   
system = init.read_xml(filename=sysfile)
r=13.0
"""
table=pair.table(width=1000)
r_cut=5*distance_scale
table.pair_coeff.set('V', 'V', func=dirac_potential, rmax=r+1, rmin = r-1, coeff=dict(r0=r))
table.pair_coeff.set('CH2', 'CH3', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict( sigma = distance_scale*3.86, epsilon = 78*kb))
table.pair_coeff.set('CH2', 'CH2', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict(sigma =distance_scale* 3.96, epsilon = 56*kb ))
table.pair_coeff.set('CH2', 'S', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict(sigma = distance_scale*4.21, epsilon =84*kb ))
table.pair_coeff.set('CH2', 'Au', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict(sigma = distance_scale*3.54, epsilon =88*kb ))
table.pair_coeff.set('CH2', 'V', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict(sigma = distance_scale*3.54, epsilon =88*kb ))
table.pair_coeff.set('CH3', 'CH3', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict( sigma = distance_scale*3.76, epsilon = 108*kb))
table.pair_coeff.set('S', 'CH3', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict( sigma =distance_scale* 4.11, epsilon = 117*kb))
table.pair_coeff.set('Au', 'CH3', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict(sigma = distance_scale* 3.54, epsilon = 108*kb ))
table.pair_coeff.set('V', 'CH3', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict(sigma = distance_scale* 3.54, epsilon = 108*kb ))
table.pair_coeff.set('S', 'S', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict( sigma =distance_scale* 4.45,epsilon = 126*kb))
table.pair_coeff.set('S', 'Au', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict( sigma =  distance_scale*2.65,epsilon = 2795*kb))
table.pair_coeff.set('S', 'V', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict( sigma =  distance_scale*2.65,epsilon = 2795*kb))
table.pair_coeff.set('Au', 'Au', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict(sigma = 1* distance_scale,epsilon = 696.379))
table.pair_coeff.set('Au', 'V', func=lj, rmax=r_cut, rmin = 0.001, coeff=dict( sigma = 1* distance_scale,epsilon = 696.379))
"""

lj=pair.lj(r_cut=2.5 * distance_scale)
lj.pair_coeff.set('CH2','CH3', sigma = distance_scale*3.86, epsilon = 78*kb, alpha=1) 
lj.pair_coeff.set('CH2','CH2',sigma =distance_scale* 3.96, epsilon = 56*kb, alpha = 1)
lj.pair_coeff.set('CH2','S', sigma = distance_scale*4.21, epsilon =84*kb, alpha=1)
lj.pair_coeff.set('CH2','Au',sigma = distance_scale*3.54, epsilon =88*kb, alpha = 1)
lj.pair_coeff.set('CH2','V',sigma = distance_scale*3.54, epsilon =88*kb, alpha = 1)
lj.pair_coeff.set('CH3','CH3', sigma = distance_scale*3.76, epsilon = 108*kb, alpha = 1)
lj.pair_coeff.set('CH3','S', sigma =distance_scale* 4.11, epsilon = 117*kb, alpha=1)
lj.pair_coeff.set('CH3','Au', sigma = distance_scale* 3.54, epsilon = 108*kb, alpha=1)
lj.pair_coeff.set('CH3','V', sigma =  distance_scale*3.54,epsilon = 108*kb, alpha=1)
lj.pair_coeff.set('S','S',  sigma =distance_scale* 4.45,epsilon = 126*kb, alpha=1)
lj.pair_coeff.set('S','Au', sigma =  distance_scale*2.65,epsilon = 2795*kb,alpha=1)
lj.pair_coeff.set('S','V',  sigma = distance_scale*2.65,epsilon = 2795*kb, alpha=1)
lj.pair_coeff.set('Au', 'Au', sigma = 1*distance_scale,epsilon = 696.379, alpha=0, r_cut=distance_scale)
lj.pair_coeff.set('V', 'Au', sigma = 1* distance_scale,epsilon = 696.379, alpha=0, r_cut=distance_scale)
lj.pair_coeff.set('V', 'V', sigma = 1* distance_scale,epsilon = 696.379, alpha=0, r_cut=distance_scale)


harmonica = angle.harmonic()
#harmonica.set_coeff('Au-S-CH2',k=(62.5e3)*kb, t0=1.7453)
harmonica.set_coeff('S-CH2-CH2',k=(62.5e3)*kb, t0=1.996)
harmonica.set_coeff('CH2-CH2-CH2',k=(62.5e3)*kb, t0=1.911)
harmonica.set_coeff('CH2-CH2-CH3',k=(62.5e3)*kb, t0=1.911)

harmonicad1 = dihedral.harmonic()
harmonicad2 = dihedral.harmonic()
harmonicad3 = dihedral.harmonic()
harmonicad4 = dihedral.harmonic()
harmonicad5 = dihedral.harmonic()
harmonicad1.set_coeff('phi1',k=(2.363e3)*kb, d=1, n=1)
harmonicad2.set_coeff('phi1',k=(1.578e3)*kb, d=1, n=2)
harmonicad3.set_coeff('phi1',k=(2.55e3)*kb, d=1, n=3)
harmonicad4.set_coeff('phi1',k=(.789e3)*kb , d=1, n=4)
harmonicad5.set_coeff('phi1',k=(.4735e3)*kb, d=1, n=5)


harmonic = bond.harmonic()
harmonic.set_coeff('S-C', k=10.0, r0=distance_scale* 1.82)
harmonic.set_coeff('C-C', k=10.0, r0=distance_scale*1.53)
harmonic.set_coeff('V-V', k=1000.0, r0=r)
#harmonic.set_coeff('V-V', k=10000.0, r0=13)
#first = group.tags(name="first",tag_min=6)
#second = group.tags(name="second",tag_min=37)
"""
table=bond.table(width=1000)
table.bond_coeff.set('V-V',func=dirac_bond,rmin=0, rmax=100,coeff=dict(kappa=300, r0=r))
table.bond_coeff.set('C-C',func=harmonic,rmin=0, rmax=100,coeff=dict(kappa=10.0,r0=distance_scale*1.53))
table.bond_coeff.set('S-C',func=harmonic,rmin=0, rmax=100,coeff=dict(kappa=10.0,r0=distance_scale*1.82))
"""
t=0
nonrigid=group.nonrigid()
rigid = group.rigid()
nlist.reset_exclusions(exclusions=['body','bond','angle','dihedral'])
nlist.set_params(check_period=1)
fire= integrate.mode_minimize_fire(group=nonrigid, dt=0.00005)
fire_r = integrate.mode_minimize_rigid_fire(group=rigid,dt=0.00005)
run(1e5)
t+=1e5
fire = integrate.mode_minimize_fire(group=nonrigid,dt=0.0005)
fire_r = integrate.mode_minimize_rigid_fire(group=rigid,dt=0.0005)
run(1e5)
t+=2e6
fire = integrate.mode_minimize_fire(group=nonrigid,dt=0.005)
fire_r = integrate.mode_minimize_rigid_fire(group=rigid,dt=0.005)
run(1e5)
t+=2e6
del fire
del fire_r

run_time=1e5
dumps_per_distance=1e4
n=dumps_per_distance
dump_time=run_time/n
xml = dump.xml(filename="atoms.dump", period=dump_time, restart=True, time_step=False)
xml.set_params(all=True)
integrate.mode_standard(dt=0.001)
nonrigid_integrator=integrate.nvt(group=nonrigid,T=real_temp*kb,tau=.65)
rigid_integrator=integrate.nvt_rigid(group=rigid,T=real_temp*kb,tau=.65)
dcd = dump.dcd(filename='mixture.dcd', period=1000)
reals=np.array([13.0])
rs=np.array([13.0])
H=np.array([-3.0]) 
while(r>6.0):
#while(r>12.9):
#table.pair_coeff.set('V', 'V', func=dirac_potential, rmax=29, rmin = 0, coeff=dict(r0=r))
    harmonic.set_coeff('V-V', k=50.0, r0=r)
    integrate.mode_standard(dt=.001)
    rigid_integrator.enable()
    nonrigid_integrator.enable()
    #table.bond_coeff.set('V-V',func=dirac_bond, rmin=0, rmax=100,coeff=dict(kappa=300,r0=r))            
    for c in range(int(n)):
        run(dump_time)
        t+=dump_time
        z=vd.v_pos_matrix('atoms.dump')
        w=np.subtract(z[2],z[1])
        w=u.length(w)
    	H=np.append(H,[w],axis=0)
    r-=(.2*distance_scale)
    harmonic.set_coeff('V-V', k=10000.0, r0=r)
    rs=np.append(rs,[r],axis=0)
    run(4e4)
    z=vd.v_pos_matrix('atoms.dump')
    w=np.subtract(z[2],z[1])
    w=u.length(w)
    reals=np.append(reals,[w],axis=0)
    rigid_integrator.disable()
    nonrigid_integrator.disable()
    fire = integrate.mode_minimize_fire(group=nonrigid,dt=0.005)
    fire_r = integrate.mode_minimize_rigid_fire(group=rigid,dt=0.005)
    run(1e4)
    del fire
    del fire_r
    z=vd.v_pos_matrix('atoms.dump')
    w=np.subtract(z[2],z[1])
    w=u.length(w)
    reals=np.append(reals,[w],axis=0)
with open('H.txt','a') as myfile:
	for i in range(1,len(H)):
		myfile.write(str(H[i])+'\n')
print rs
print '\n'
print reals
