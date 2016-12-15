"""
Created on Tue Jul  5 10:44:23 2016
@author: waltmann
"""

from __future__ import division
#from hoomd_script import *
import math
import sys
sys.path.insert(0, "../analysis/hydrocarbon/")
print sys.path
print '\n\n\n\n'
import hoomd
from hoomd import deprecated
from hoomd import md
#import util as u
import hoomdinit as hi
import numpy as np
#import moment_func as mf
kb=8.61733e-5*30
distance_scale=1.0/3.96
real_temp=330
bond_strength=30
#8.18+1.97+1.25
npRadius=1.14
repeats=12
l=.128*repeats+.2
lamda=l/npRadius
opmtau=(3*lamda+1)**(1/3)
ocmtau=-1*(1+lamda)/2+(((1+lamda)/2)**2+(6*lamda+2)/(1+lamda))**(.5)


hoomd.context.initialize("--mode=cpu")
sysfile= "Au140HC_r20.xml" 
system = hoomd.deprecated.init.read_xml(filename=sysfile)
		
r=20
steps=10
instep=(r-(ocmtau-4)*distance_scale)/steps


rigid = hoomd.md.constrain.rigid()

hi.set_up_rigid(system,sysfile,rigid,distance_scale=1/distance_scale,)

rigid.create_bodies(False)

nl = hoomd.md.nlist.cell()
nl.set_params(check_period=1)
lj=hoomd.md.pair.lj(r_cut=3/ distance_scale,nlist=nl)
nl.reset_exclusions(exclusions=['body','bond','angle','dihedral'])
lj.pair_coeff.set('CH2','CH3', sigma = distance_scale*3.86, epsilon = 78*kb, alpha=1) 
lj.pair_coeff.set('CH2','CH2',sigma =distance_scale* 3.96, epsilon = 56*kb, alpha = 1)
lj.pair_coeff.set('CH2','S', sigma = distance_scale*4.21, epsilon =84*kb, alpha=1)
lj.pair_coeff.set('CH2','Au',sigma = distance_scale*3.54, epsilon =88*kb, alpha = 1)
lj.pair_coeff.set('CH2','V1',sigma = distance_scale*3.54, epsilon =88*kb, alpha = 1)
lj.pair_coeff.set('CH3','CH3', sigma = distance_scale*3.76, epsilon = 108*kb, alpha = 1)
lj.pair_coeff.set('CH3','S', sigma =distance_scale* 4.11, epsilon = 117*kb, alpha=1)
lj.pair_coeff.set('CH3','Au', sigma = distance_scale* 3.54, epsilon = 108*kb, alpha=1)
lj.pair_coeff.set('CH3','V1', sigma =  distance_scale*3.54,epsilon = 108*kb, alpha=1)
lj.pair_coeff.set('S','S',  sigma =distance_scale* 4.45,epsilon = 126*kb, alpha=1)
lj.pair_coeff.set('S','Au', sigma =  distance_scale*2.65,epsilon = 2795*kb,alpha=1)
lj.pair_coeff.set('S','V1',  sigma = distance_scale*2.65,epsilon = 2795*kb, alpha=1)
lj.pair_coeff.set('Au', 'Au', sigma = 1*distance_scale,epsilon = 696.379, alpha=0, r_cut=distance_scale)
lj.pair_coeff.set('V1', 'Au', sigma = 1* distance_scale,epsilon = 696.379, alpha=0, r_cut=distance_scale)
lj.pair_coeff.set('V1', 'V1', sigma = 1* distance_scale,epsilon = 696.379, alpha=0, r_cut=distance_scale)
lj.pair_coeff.set('V2', 'Au', sigma = 1* distance_scale,epsilon = 696.379, alpha=0, r_cut=distance_scale)
lj.pair_coeff.set('V2', 'V2', sigma = 1* distance_scale,epsilon = 696.379, alpha=0, r_cut=distance_scale)
lj.pair_coeff.set('S','V2',  sigma = distance_scale*2.65,epsilon = 2795*kb, alpha=1)
lj.pair_coeff.set('CH3','V2', sigma =  distance_scale*3.54,epsilon = 108*kb, alpha=1)
lj.pair_coeff.set('CH2','V2',sigma = distance_scale*3.54, epsilon =88*kb, alpha = 1)
lj.pair_coeff.set('V2', 'V1', sigma = 1* distance_scale,epsilon = 696.379, alpha=0, r_cut=distance_scale)

harmonica = hoomd.md.angle.harmonic()
#harmonica.set_coeff('Au-S-CH2',k=(62.5e3)*kb, t0=1.7453)
harmonica.angle_coeff.set('S-CH2-CH2',k=(62.5e3)*kb, t0=1.996)
harmonica.angle_coeff.set('CH2-CH2-CH2',k=(62.5e3)*kb, t0=1.911)
harmonica.angle_coeff.set('CH2-CH2-CH3',k=(62.5e3)*kb, t0=1.911)

harmonicad1 = hoomd.md.dihedral.harmonic()
harmonicad2 = hoomd.md.dihedral.harmonic()
harmonicad3 = hoomd.md.dihedral.harmonic()
harmonicad4 = hoomd.md.dihedral.harmonic()
harmonicad5 = hoomd.md.dihedral.harmonic()
harmonicad6 = hoomd.md.dihedral.harmonic()
harmonicad1.dihedral_coeff.set('phi1',k=(2.363e3)*kb, d=1, n=1)
harmonicad2.dihedral_coeff.set('phi1',k=(1.578e3)*kb, d=1, n=2)
harmonicad3.dihedral_coeff.set('phi1',k=(2.55e3)*kb, d=1, n=3)
harmonicad4.dihedral_coeff.set('phi1',k=(.789e3)*kb , d=1, n=4)
harmonicad5.dihedral_coeff.set('phi1',k=(.4735e3)*kb, d=1, n=5)
harmonicad6.dihedral_coeff.set('phi1',k=-6.119,d=1,n=0)

harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('S-C', k=12000.0, r0=distance_scale* 1.82)
harmonic.bond_coeff.set('C-C', k=12000.0, r0=distance_scale*1.53)
harmonic.bond_coeff.set('V-V', k=bond_strength, r0=r)


t=0
nonrigid=hoomd.group.nonrigid()
centers = hoomd.group.rigid_center()
particles= hoomd.group.union(name='paricles',a=centers,b=nonrigid) 


#fire= md.integrate.mode_minimize_fire(group=particles, dt=0.001)
#hoomd.run(1e5)
"""
t+=1e5
fire = hoomd.md.integrate.mode_minimize_fire(group=particles,dt=0.0005)
hoomd.run(1e5)
t+=2e6
fire = hoomd.md.integrate.mode_minimize_fire(group=particles,dt=0.005)
hoomd.run(1e5)
t+=2e6
del fire
"""




run_time=1e3
dumps_per_distance=1e2
n=dumps_per_distance
logger = hoomd.analyze.log(filename='mylog.log', period=run_time/n, quantities=['temperature','potential_energy','kinetic_energy','volume','pressure','pair_lj_energy','bond_harmonic_energy','angle_harmonic_energy'])
dump_time=run_time/n
xml = hoomd.deprecated.dump.xml(filename="atoms.dump", period=dump_time, group= hoomd.group.all(),time_step=False, restart=True)
xml.set_params(all=True)
hoomd.md.integrate.mode_standard(dt=0.001)
nonrigid_integrator=hoomd.md.integrate.nvt(group=particles,kT=real_temp*kb,tau=.65)
dcd = hoomd.dump.dcd(filename='mixture.dcd', period=1000, overwrite=True)
reals=np.array([r])
rs=np.array([r])
H=np.array([-3.0]) 

while(r>(ocmtau-4)*distance_scale):
    hoomd.md.integrate.mode_standard(dt=.001)
    logger.enable()           
    for c in range(int(n)):
        	hoomd.run(dump_time)
        	t+=dump_time
        	z=hi.v_pos_matrix('atoms.dump')
        	w=hi.part_distance(z[2],z[1],'atoms.dump')
    		H=np.append(H,[w],axis=0)

    r-=instep
    harmonic.bond_coeff.set('V-V', k=bond_strength, r0=r)
    rs=np.append(rs,[r],axis=0)
    logger.disable()
    z=hi.v_pos_matrix('atoms.dump')
    w=hi.part_distance(z[2],z[1],'atoms.dump')
    reals=np.append(reals,[w],axis=0)
    H=np.append(H,[-3.0],axis=0)
    harmonic.bond_coeff.set('V-V', k=bond_strength, r0=r)
    try:    
        hoomd.run(2e4)
    except(RuntimeError):
        print("Failed at " + str(r))
        break
    z=hi.v_pos_matrix('atoms.dump')
    w=hi.part_distance(z[2],z[1],'atoms.dump')
    reals=np.append(reals,[w],axis=0)
H=H[1:]
with open('H.txt') as myfile:
	for i in range(1,len(H)):
		myfile.write(str(H[i])+'\n')
with open('R.txt') as rfile:
    for i in range(0,len(rs)):
        rfile.write(str(rs[i]) + '\n')
print 'where its supposed to be'
print rs
print '\n'
print 'where it is'
print reals
