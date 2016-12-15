# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 09:45:24 2016

@author: waltmann
"""

import hoomd
from hoomd import deprecated
from hoomd import md
import numpy as np
import hoomdinit as hi

kb= .002585199 ##eV/K (8.61733e-5 * 30)
distance_scale=1.0/3.96
hoomd.context.initialize("--mode=cpu")
sysfile='Au201HC.xml'
system=hoomd.deprecated.init.read_xml(sysfile)

"""
for p in system.particles:
	if (p.type[:1]=='V' or p.type=='Au' or p.type=='S'):
		p.mass=1
"""
## energies*30 distances/3.96 temp=375*30
kb30=0.002585
dist_scale=1/3.96
run_time=1000000
traj_period=1000
log_period=1000
dump_period=10000

harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('C-C', k=12000, r0=dist_scale*1.53)
harmonic.bond_coeff.set('S-C', k=12000, r0=dist_scale*1.82)

harmonica = hoomd.md.angle.harmonic()
harmonica.angle_coeff.set('CH2-CH2-CH2',k=62.5e3*kb,t0=1.911)
harmonica.angle_coeff.set('CH2-CH2-CH3',k=62.5e3*kb,t0=1.911)
harmonica.angle_coeff.set('S-CH2-CH2',k=62.5e3*kb,t0=1.996)


harmonicad1 = hoomd.md.dihedral.harmonic()
harmonicad2 = hoomd.md.dihedral.harmonic()
harmonicad3 = hoomd.md.dihedral.harmonic()
harmonicad1.dihedral_coeff.set('phi1',k=(1.835604), d=1, n=1)
harmonicad2.dihedral_coeff.set('phi1',k=(.352504), d=-1, n=2)
harmonicad3.dihedral_coeff.set('phi1',k=(4.0914), d=1, n=3)

rigid = hoomd.md.constrain.rigid()

hi.set_up_rigid(system,sysfile,rigid,distance_scale=1/distance_scale,)

rigid.create_bodies(False)


nl = hoomd.md.nlist.cell()
nl.set_params(check_period=1)
nl.reset_exclusions(exclusions=['body','bond','angle','dihedral'])
lj=hoomd.md.pair.lj(r_cut=3/distance_scale,nlist=nl)
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
lj.pair_coeff.set('Au', 'Au', sigma = 1,epsilon = 1*kb, alpha=0)
lj.pair_coeff.set('V1', 'Au', sigma = 1,epsilon = 1*kb, alpha=0)
lj.pair_coeff.set('V1', 'V1', sigma = 1,epsilon = 1*kb, alpha=0)

nonrigid=hoomd.group.nonrigid()


centers=hoomd.group.rigid_center()

particles=hoomd.group.union(name='particles',a=centers,b=nonrigid)
dcd=hoomd.dump.dcd(filename='trajectory.dcd',period=traj_period,overwrite=True)
xml=hoomd.deprecated.dump.xml(group=hoomd.group.all(),filename='atoms.dump',all=True,period=dump_period)
#fire= hoomd.md.integrate.mode_minimize_fire(group=particles, dt=0.00005)
#hoomd.run(1e5)
"""
fire = hoomd.md.integrate.mode_minimize_fire(group=particles,dt=0.001)
hoomd.run(1e5)

del fire_n
"""

hoomd.md.integrate.mode_standard(dt=0.001, aniso=True)
Tvar=hoomd.variant.linear_interp(points=[(0,20*kb30),(run_time,500*kb30)])
hoomd.md.integrate.nvt(group=particles, kT=400*kb, tau=0.5)
#hoomd.md.integrate.langevin(group=particles, kT=400*kb, seed=1)
#hoomd.md.integrate.npt(group=particles, T=Tvar, tau=0.5)


logger=hoomd.analyze.log(filename='mylog.log',quantities=['temperature','potential_energy','kinetic_energy','volume','pressure'],period=log_period,overwrite=True)
#zeroer= hoomd.md.update.zero_momentum(period=10)
hoomd.run(run_time)







