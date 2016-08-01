from __future__ import division
from hoomd_script import *
import math
import sys
import os
import numpy as np
##################################################
##Fundamental constants
################################################
kb= 8.61733e-5 * 30##eV/K
c=2.99792e6  #A/ps
anorm=14.0266 ##mass normalization factor for unified atoms
amu=9.31494e8  ##eV/(c^2)
distance_scale=1/3.96 ###distance unit/angstroms 
######################################################
## Get Initial XML FILE
######################################################
system = init.read_xml(filename="Lat2by2by2_25_Fcc_minimized.xml")
######################################################
#  Parameters That Can Be Changed
######################################################
real_temp=300
timestep=0.0001  ### never greater than .001
run_time=1000000
log_period=run_time/1000
ndumpfiles=5
nframes=500
box_start=50.0
box_end=22
r_start=25.0 * np.sqrt(2)/2.0
steps=20
resize_time=(run_time/4)
######################################################
# Bond, Angle and Dihedral Setup
######################################################
harmonic = bond.harmonic()
harmonic.set_coeff('S-C', k=300.0, r0=distance_scale* 1.82)
harmonic.set_coeff('C-C', k=300.0, r0=distance_scale*1.53)
harmonic.set_coeff('V-V', k=3000.0, r0=r_start)
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
######################################################
#attraction and repulsion parameters
######################################################
#force field setup
lj=pair.lj(r_cut=distance_scale*2.5)

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
"""
lj.pair_coeff.set('S','Au', sigma =  distance_scale*2.65,epsilon = 2795*kb,alpha=1)
lj.pair_coeff.set('S','V',  sigma = distance_scale*2.65,epsilon = 2795*kb, alpha=1)
lj.pair_coeff.set('Au', 'Au', sigma = 1*distance_scale,epsilon = 696*kb, alpha=1)
lj.pair_coeff.set('V', 'Au', sigma = 1*distance_scale,epsilon = 696*kb, alpha=1)
lj.pair_coeff.set('V', 'V', sigma = 1*distance_scale,epsilon = 696*kb, alpha=1)
"""
lj.pair_coeff.set('S','Au', sigma =  distance_scale*2.65,epsilon = 0,alpha=1)
lj.pair_coeff.set('S','V',  sigma = distance_scale*2.65,epsilon = 0, alpha=1)
lj.pair_coeff.set('Au', 'Au', sigma = 1*distance_scale,epsilon = 0, alpha=1)
lj.pair_coeff.set('V', 'Au', sigma = 1*distance_scale,epsilon = 0, alpha=1)
lj.pair_coeff.set('V', 'V', sigma = 1*distance_scale,epsilon = 0, alpha=1)
######################################################
# Group nonrigid particles and minimize nonrigid energy
######################################################
rigid = group.rigid()
nonrigid = group.nonrigid()
nlist.reset_exclusions(exclusions=['body','bond','angle','dihedral'])
nlist.set_params(check_period=1)
"""
fire= integrate.mode_minimize_rigid_fire(group=rigid, dt=0.00005)
fire_r = integrate.mode_minimize_fire(group=nonrigid,dt=0.00005)
run(1e5)
fire = integrate.mode_minimize_fire(group=nonrigid,dt=0.0005)
fire_r = integrate.mode_minimize_rigid_fire(group=rigid, dt=0.0005)
run(1e5)
fire = integrate.mode_minimize_fire(group=nonrigid,dt=0.005)
fire_r = integrate.mode_minimize_rigid_fire(group=rigid, dt=0.005)
run(1e5)
del fire
del fire_r
"""
#############################################
#Zero The momentum
#############################################
zeroer= update.zero_momentum(period=500000)
#change the neighborlist check time for boost in performance
nlist.set_params(check_period=1)
#set integrate back to standard dt
####################################################
integrate.mode_standard(timestep)
####################################################
#       Log Kinetic Energy and Temperature of System
####################################################
logger = analyze.log(filename='mylog.log', period=log_period, quantities=['temperature','potential_energy','kinetic_energy','volume','pressure'])
pressure = analyze.log(filename='pressure.log', period=log_period, quantities=['pressure','pressure_xx','pressure_yy','pressure_zz','pressure_xy','pressure_yz','pressure_xz'])                            
######################################################
#       Dump Files
######################################################
#xml dump files
xml = dump.xml(filename="atoms.dump", period=run_time/ndumpfiles)
xml.set_params(all=True)
# dump a .dcd file for the trajectory
dcd = dump.dcd(filename='mixture.dcd', period=run_time/nframes)
######################################################
# make group of rigid particles and set up integrators
######################################################
rigid_integrator=integrate.nvt_rigid(group=rigid, T=real_temp*kb, tau=0.65)
#rigid_integrator=integrate.nvt_rigid(group=rigid, T=variant.linear_interp([(0,20*kb),(run_time,500*kb)]), tau=0.65)
nonrigid_integrator=integrate.nvt(group=nonrigid, T=real_temp*kb, tau=0.65)
box_resize = update.box_resize(L=variant.linear_interp([(0,box_start),(resize_time,box_end)]))
                              

####################################################
######################################################
#   Main Run
######################################################
time=0  
r=r_start
while(time<resize_time):
    r-=(r_start-((box_end/box_start)*r_start))/steps
    harmonic.set_coeff('V-V', k=30.0, r0=r)
    run(resize_time/steps)
    print '\n'
    print r 
    print '\n'
    time+=resize_time/steps

harmonic.set_coeff('V-V', k=3000.0, r0=r)

run(run_time)
