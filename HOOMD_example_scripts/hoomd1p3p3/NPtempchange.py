from __future__ import division
from hoomd_script import *
import math
import sys
import os
##################################################
##Fundamental constants
################################################
kb= .002585199 ##eV/K (8.61733e-5 * 30)
c=2.99792e6  #A/ps
anorm=14.0266 ##mass normalization factor for unified atoms
amu=9.31494e8  ##eV/(c^2)
distance_scale=1/3.96 ###distance unit/angstroms 
######################################################
## Get Initial XML FILE
######################################################
system = init.read_xml(filename="Au201HC.xml")
######################################################
#  Parameters That Can Be Changed
######################################################
low_temp=250
real_temp=280
real_temp_high=350
timestep=0.001  ### never greater than .001
run_time=500000
log_period=10000
######################################################
# Bond, Angle and Dihedral Setup
######################################################
harmonic = bond.harmonic()
harmonic.set_coeff('S-C', k=12000.0, r0=distance_scale* 1.82)
harmonic.set_coeff('C-C', k=12000.0, r0=distance_scale*1.53)

harmonica = angle.harmonic()
#harmonica.set_coeff('Au-S-CH2',k=(62.5e3)*kb, t0=1.7453)
harmonica.set_coeff('S-CH2-CH2',k=(62.5e3)*kb, t0=1.996)
harmonica.set_coeff('CH2-CH2-CH2',k=(62.5e3)*kb, t0=1.911)
harmonica.set_coeff('CH2-CH2-CH3',k=(62.5e3)*kb, t0=1.911)


harmonicad1 = dihedral.harmonic()
harmonicad2 = dihedral.harmonic()
harmonicad3 = dihedral.harmonic()
harmonicad1.set_coeff('phi1',k=(1.835604), d=1, n=1)
harmonicad2.set_coeff('phi1',k=(.352504), d=-1, n=2)
harmonicad3.set_coeff('phi1',k=(4.0914), d=1, n=3)
######################################################
#attraction and repulsion parameters
######################################################
#force field setup
lj=pair.lj(r_cut=3/distance_scale)

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
rigid_integrator=integrate.nvt_rigid(group=rigid, T=500*kb, tau=0.65)
nonrigid_integrator=integrate.nvt(group=nonrigid, T=500*kb, tau=0.65)
######################################################
# Equilibrate System and shrink box
######################################################
######################################################
integrate.mode_standard(dt=0.0005)
#box_resize = update.box_resize(Lx=variant.linear_interp([(0,50),(5e5,55)]),Ly=variant.linear_interp([(0,50),(5e5,55)]),Lz=variant.linear_interp([(0,50),(5e5,55)]))
#run(10000)
#############################################
#Zero The momentum
#############################################
#change the neighborlist check time for boost in performance
nlist.set_params(check_period=3)
####################################################
#set integrate back to standard dt
####################################################
integrate.mode_standard(timestep)
####################################################
#       Log Kinetic Energy and Temperature of System
####################################################
logger = analyze.log(filename='mylog.log', period=log_period, quantities=['temperature','potential_energy','kinetic_energy','volume','pressure','pair_lj_energy','bond_harmonic_energy','angle_harmonic_energy'])
pressure = analyze.log(filename='pressure.log', period=log_period, quantities=['pressure','pressure_xx','pressure_yy','pressure_zz','pressure_xy','pressure_yz','pressure_xz'])                            
######################################################
#       Dump Files
######################################################
#xml dump files
xml = dump.xml(filename="atoms.dump", period=run_time/100)
xml.set_params(all=True)
# dump a .dcd file for the trajectory
dcd = dump.dcd(filename='mixture.dcd', period=5000)
######################################################
#   Main Run
######################################################  
tvar=variant.linear_interp([(0,500*kb),(run_time,500*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,500*kb),(run_time,450*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,450*kb),(run_time,450*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,450*kb),(run_time,400*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,400*kb),(run_time,400*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,400*kb),(run_time,350*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,350*kb),(run_time,350*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,350*kb),(run_time,330*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,330*kb),(run_time,330*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,330*kb),(run_time,320*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,320*kb),(run_time,320*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,320*kb),(run_time,310*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,310*kb),(run_time,310*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,310*kb),(run_time,300*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,300*kb),(run_time,300*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,300*kb),(run_time,290*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,290*kb),(run_time,290*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,290*kb),(run_time,250*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,250*kb),(run_time,250*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,250*kb),(run_time,200*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,200*kb),(run_time,200*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,200*kb),(run_time,150*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,150*kb),(run_time,150*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,150*kb),(run_time,100*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,100*kb),(run_time,100*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,100*kb),(run_time,50*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,50*kb),(run_time,50*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,50*kb),(run_time,1*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
tvar=variant.linear_interp([(0,1*kb),(run_time,1*kb)])
rigid_integrator.set_params(T=tvar)
nonrigid_integrator.set_params(T=tvar)
run(run_time)
