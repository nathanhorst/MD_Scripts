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
distance_scale=1/3.96 ###distance unit/angstroms 
######################################################
## Get Initial XML FILE
######################################################
system = init.read_xml(filename="scaled_np.xml")
  
######################################################
#  Parameters That Can Be Changed
######################################################
  
real_temp=100
timestep=0.001  ### never greater than .001
run_time=1000000
log_period=10000
  
######################################################
# Bond, Angle and Dihedral Setup
######################################################
 
harmonic = bond.harmonic()
harmonic.set_coeff('S-C', k=10.0, r0=distance_scale* 1.82)
harmonic.set_coeff('C-C', k=10.0, r0=distance_scale*1.53)
 
harmonica = angle.harmonic()
harmonica.set_coeff('Au-S-CH2',k=(62.5e3)*kb, t0=1.7453)
harmonica.set_coeff('S-CH2-CH2',k=(62.5e3)*kb, t0=1.996)
harmonica.set_coeff('CH2-CH2-CH2',k=(62.5e3)*kb, t0=1.911)
harmonica.set_coeff('CH2-CH2-CH3',k=(62.5e3)*kb, t0=1.911)

"""
def torsion(theta, k, n,d):
    V = 0.5 * k * (1+d(cos(n*theta)));
    F = -0.5*k*d*n * sin(n *theta);
    return (V, F)
     
dtable = dihedral.table(width=1000)
dtable.dihedral_coeff.set('phi1', func=torsion, coeff=dict(k=(2.363e3)*kb, d=-1, n=1))
dtable.dihedral_coeff.set('phi2', func=torsion, coeff=dict(k=(1.578e3)*kb, d=1, n=2))
dtable.dihedral_coeff.set('phi3', func=torsion, coeff=dict(k=(2.55e3)*kb, d=-1, n=3))
dtable.dihedral_coeff.set('phi4', func=torsion, coeff=dict(k=(.789e3)*kb , d=1, n=4))
dtable.dihedral_coeff.set('phi5', func=torsion, coeff=dict(k=(.4735e3)*kb, d=-1, n=5))
"""

harmonicad = dihedral.harmonic()
harmonicad.set_coeff('phi1',k=(2.363e3)*kb, d=1, n=1)
harmonicad.set_coeff('phi2',k=(1.578e3)*kb, d=-1, n=2)
harmonicad.set_coeff('phi3',k=(2.55e3)*kb, d=1, n=3)
harmonicad.set_coeff('phi4',k=(.789e3)*kb , d=-1, n=4)
harmonicad.set_coeff('phi5',k=(.4735e3)*kb, d=1, n=5)

######################################################
#attraction and repulsion parameters
######################################################
 
#force field setup
lj=pair.lj(r_cut=distance_scale*10)

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
lj.pair_coeff.set('Au', 'Au', sigma = 1,epsilon = 1*kb, alpha=0)
lj.pair_coeff.set('V', 'Au', sigma = 1,epsilon = 1*kb, alpha=0)
lj.pair_coeff.set('V', 'V', sigma = 1,epsilon = 1*kb, alpha=0)

######################################################
# Group nonrigid particles and minimize nonrigid energy
######################################################
rigid = group.rigid()
nonrigid = group.nonrigid()
nlist.reset_exclusions(exclusions=['body','bond','angle','dihedral'])
nlist.set_params(check_period=1)
fire= integrate.mode_minimize_rigid_fire(group=rigid, dt=0.00005)
fire_r = integrate.mode_minimize_fire(group=nonrigid,dt=0.00005)
run(5e5)
fire = integrate.mode_minimize_fire(group=nonrigid,dt=0.0005)
fire_r = integrate.mode_minimize_rigid_fire(group=rigid, dt=0.0005)
run(2e6)
fire = integrate.mode_minimize_fire(group=nonrigid,dt=0.005)
fire_r = integrate.mode_minimize_rigid_fire(group=rigid, dt=0.005)
run(2e6)
del fire
del fire_r
######################################################
# make group of rigid particles and set up integrators
######################################################
  
#integrate of rigid and nonrigid
 
rigid_integrator=integrate.nvt_rigid(group=rigid, T=variant.linear_interp([(0,20*kb),(run_time,500*kb)]), tau=0.65)
nonrigid_integrator=integrate.nvt(group=nonrigid, T=variant.linear_interp([(0,20*kb),(run_time,500*kb)]), tau=0.65)
  
######################################################
# Equilibrate System and shrink box
######################################################
######################################################
#set integrate low so that system can equilibrate
#set the check period very low to let the system equilibrate
 
integrate.mode_standard(dt=0.0005)

#box_resize = update.box_resize(Lx=variant.linear_interp([(0,50),(5e5,55)]),Ly=variant.linear_interp([(0,50),(5e5,55)]),Lz=variant.linear_interp([(0,50),(5e5,55)]))
#run(10000)
#############################################
#Zero The momentum
#############################################
  
zeroer= update.zero_momentum(period=500000)
#change the neighborlist check time for boost in performance
nlist.set_params(check_period=3)
  
####################################################
#set integrate back to standard dt
####################################################
integrate.mode_standard(timestep)
  
####################################################
#       Log Kinetic Energy and Temperature of System
####################################################
  
logger = analyze.log(filename='mylog.log', period=log_period, quantities=['temperature','potential_energy','kinetic_energy','volume','pressure'])
#pressure = analyze.log(filename='pressure.log', period=log_period, quantities=['pressure',

 #                       'pressure_xx','pressure_yy','pressure_zz','pressure_xy',
  #                      'pressure_yz','pressure_xz'])                            
######################################################
#       Dump Files
######################################################
  
#xml dump files
xml = dump.xml(filename="atoms.dump", period=run_time/96)
xml.set_params(all=True)
# dump a .dcd file for the trajectory
dcd = dump.dcd(filename='mixture.dcd', period=5000)
######################################################
#   Main Run
######################################################  
run(run_time)
