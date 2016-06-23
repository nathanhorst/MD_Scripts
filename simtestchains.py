# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 15:48:25 2016

@author: waltmann
"""

from __future__ import division
from hoomd_script import *
import math
import sys
import os

 
 
##################################################
##Fundamental constants
################################################
 
kb= 8.61733e-5 ##eV/K
c=2.99792e6  #A/ps
anorm=14.0266 ##mass normalization factor for unified atoms
amu=9.31494e8  ##eV/(c^2)
distance_scale=1 ###distance unit/angstroms 
######################################################
## Get Initial XML FILE
######################################################
system = init.read_xml(filename="trans_chain.xml")
  
######################################################
#  Parameters That Can Be Changed
######################################################
  
real_temp=1000
timestep=0.0001  ### never greater than .001
run_time=5000000
log_period=2000
  
######################################################
# Bond, Angle and Dihedral Setup
######################################################
 
harmonic = bond.harmonic()
harmonic.set_coeff('S-C', k=10.0, r0=distance_scale* 1.82)
harmonic.set_coeff('C-C', k=10.0, r0=distance_scale*1.53)
"""
harmonic = bond.harmonic()
harmonic.set_coeff('S-C', k=330.0, r0=0.82)
harmonic.set_coeff('C-C', k=330.0, r0=0.84)
 """
harmonica = angle.harmonic()
#harmonica.set_coeff('Au-S-CH2',k=(62.5e3)*kb, t0=1.7453)
harmonica.set_coeff('S-CH2-CH2',k=(62.5e3)*kb, t0=1.996)
harmonica.set_coeff('CH2-CH2-CH2',k=(62.5e3)*kb, t0=1.911)
harmonica.set_coeff('CH2-CH2-CH3',k=(62.5e3)*kb, t0=1.911)
"""
harmonica = angle.harmonic()
harmonica.set_coeff('S-CH2-CH2',k=100.0, t0=1.0)
harmonica.set_coeff('CH2-CH2-CH2',k=3.0, t0=0.7851)
harmonica.set_coeff('CH2-CH2-CH3',k=3.0, t0=0.7851)
"""
"""
def torsion(theta, k, n,d):
    V = 0.5 * k * (1+d(cos(n*theta)));
    F = -0.5*k*d*n * sin(n *theta);
    return (V, F)
 """  
"""  
dtable = dihedral.table(width=1000)
dtable.dihedral_coeff.set('phi1', func=torsion, coeff=dict(k=(2.363e3)*kb, d=-1, n=1))
"""
#dtable.dihedral_coeff.set('phi2', func=torsion, coeff=dict(k=(1.578e3)*kb, d=1, n=2))
#dtable.dihedral_coeff.set('phi3', func=torsion, coeff=dict(k=(2.55e3)*kb, d=-1, n=3))
#dtable.dihedral_coeff.set('phi4', func=torsion, coeff=dict(k=(.789e3)*kb , d=1, n=4))
#dtable.dihedral_coeff.set('phi5', func=torsion, coeff=dict(k=(.4735e3)*kb, d=-1, n=5))
"""
def harmonic(theta, kappa, theta0):
    V = 0.5 * kappa * (theta-theta0)**2;
    F = -kappa*(theta-theta0);
    return (V, F)
dtable = dihedral.table(width=1000)
dtable.dihedral_coeff.set('phi1', func=harmonic, coeff=dict(kappa=330, theta_0=0.0))

harmonicad = dihedral.harmonic()
harmonicad.set_coeff('phi1',k=(2.363e3)*kb, d=-1, n=1)
harmonicad.set_coeff('phi2',k=(1.578e3)*kb, d=1, n=2)
harmonicad.set_coeff('phi3',k=(2.55e3)*kb, d=-1, n=3)
harmonicad.set_coeff('phi4',k=(.789e3)*kb , d=1, n=4)
harmonicad.set_coeff('phi5',k=(.4735e3)*kb, d=-1, n=5)
"""
harmonicad=dihedral.harmonic()
harmonicad.set_coeff('phi1',k=(2.363e3)*kb, d=1, n=1)
harmonicad.set_coeff('phi2',k=(1.578e3)*kb, d=-1, n=2)
harmonicad.set_coeff('phi3',k=(2.55e3)*kb, d=1, n=3)
harmonicad.set_coeff('phi4',k=(.789e3)*kb , d=-1, n=4)
harmonicad.set_coeff('phi5',k=(.4735e3)*kb, d=1, n=5)

######################################################
#attraction and repulsion parameters
######################################################
 
#force field setup
lj=pair.lj(r_cut=distance_scale*5)

lj.pair_coeff.set('CH2','CH3', sigma = distance_scale*3.86, epsilon = 78*kb, alpha=1) 
lj.pair_coeff.set('CH2','CH2',sigma =distance_scale* 3.96, epsilon = 56*kb, alpha = 1)
lj.pair_coeff.set('CH2','S', sigma = distance_scale*4.21, epsilon =84*kb, alpha=1)
lj.pair_coeff.set('CH3','CH3', sigma = distance_scale*3.76, epsilon = 108*kb, alpha = 1)
lj.pair_coeff.set('CH3','S', sigma =distance_scale* 4.11, epsilon = 117*kb, alpha=1)
lj.pair_coeff.set('S','S',  sigma =distance_scale* 4.45,epsilon = 126*kb, alpha=1)
#lj.pair_coeff.set('x','S',  sigma =distance_scale* 4.45,epsilon = 1260*kb, alpha=1)
#lj.pair_coeff.set('x','CH2',  sigma =distance_scale* 4.45,epsilon = 1260*kb, alpha=1)
#lj.pair_coeff.set('x','CH3',  sigma =distance_scale* 4.45,epsilon = 1260*kb, alpha=1)
#lj.pair_coeff.set('x','x',  sigma =distance_scale* 0.45,epsilon = 1206*kb, alpha=1)
"""
lj.pair_coeff.set('CH2','CH3',sigma=0.6, epsilon=1.2) 
lj.pair_coeff.set('CH2','CH2',sigma=0.6, epsilon=1.1)
lj.pair_coeff.set('CH2','S',sigma=0.6, epsilon=1.2)
lj.pair_coeff.set('CH3','CH3',sigma=0.6, epsilon=1.3)
lj.pair_coeff.set('CH3','S',sigma=0.6, epsilon=1.1)
lj.pair_coeff.set('S','S',sigma=0.6,epsilon=1.3)
"""
######################################################
# Group nonrigid particles and minimize nonrigid energy
######################################################
nonrigid = group.nonrigid()
nlist.reset_exclusions(exclusions=['body','bond','angle','dihedral'])
nlist.set_params(check_period=1)
fire_r = integrate.mode_minimize_fire(group=nonrigid,dt=0.00005)
run(5e5)
del fire_r
 #xml dump files

logger = analyze.log(filename='mylog.log', period=log_period, quantities=['temperature','potential_energy','kinetic_energy','volume','pressure'])
xml = dump.xml(filename="atoms.dump", period=run_time/96)
xml.set_params(all=True)
# dump a .dcd file for the trajectory
dcd = dump.dcd(filename='mixture.dcd', period=5000)
######################################################
#   Main Run
######################################################  
integrate.mode_standard(timestep)
nonrigid_integrator=integrate.nvt(group=nonrigid, T=variant.linear_interp([(0,20*kb ),(run_time*2 ,1000*kb)]), tau=0.65)
run(run_time*2, limit_hours=2)
