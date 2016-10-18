# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:56:34 2016

@author: nathan
"""
from __future__ import division
from hoomd_script import *

########################################
# Define constants and units of system #
# Set parameters for easily changing   #
########################################

## CONSTANTS ##
kb= 8.61733e-5 ##eV/K
kb=kb*30    ##Factors to normalize system to 1
c=2.99792e6  #A/ps
anorm=14.0266 ##mass normalization factor for unified atoms
amu=9.31494e8  ##eV/(c^2)
distance_scale=1/3.96 ###distance unit/angstroms normalization

## SYSTEM PARAMETERS ##
real_temp=330
timestep=0.001
run_time=1e5
log_period=run_time/1000
xml_period=run_time/50
traj_period=run_time/100

### Only Necessary when adjusting box size during run ###
box_start=40
box_end=40
t_start=40
t_end=40

#####################
# Initialize system #
#####################
## READ IN XML FILE ##
system = init.read_xml(filename="Au140HC.xml")

## SET HARMONIC BONDS ##
harmonic = bond.harmonic()
harmonic.set_coeff('S-C', k=12000.0, r0=distance_scale*1.82)
harmonic.set_coeff('C-C', k=12000.0, r0=distance_scale*1.53)

## SET HARMONIC ANGLES ##
harmonica = angle.harmonic()
harmonica.set_coeff('S-CH2-CH2',k=(62.5e3)*kb, t0=1.996)
harmonica.set_coeff('CH2-CH2-CH2',k=(62.5e3)*kb, t0=1.911)
harmonica.set_coeff('CH2-CH2-CH3',k=(62.5e3)*kb, t0=1.911)

## SET DIHEDRAL ANGLES ##
harmonicad = dihedral.harmonic()
harmonicad1= dihedral.harmonic()
harmonicad2= dihedral.harmonic()
harmonicad3= dihedral.harmonic()
harmonicad4= dihedral.harmonic()
harmonicad5= dihedral.harmonic()
harmonicad.set_coeff('phi1', k=-6.119, d=1, n=0)
harmonicad1.set_coeff('phi1',k=(2.363e3)*kb, d=1, n=1)
harmonicad2.set_coeff('phi1',k=(1.578e3)*kb, d=1, n=2)
harmonicad3.set_coeff('phi1',k=(2.55e3)*kb, d=1, n=3)
harmonicad4.set_coeff('phi1',k=(.789e3)*kb , d=1, n=4)
harmonicad5.set_coeff('phi1',k=(.4735e3)*kb, d=1, n=5)

## SET UP NON-BONDED INTERACTIONS ##
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
lj.pair_coeff.set('S','Au', sigma =  distance_scale*2.65,epsilon = 0,alpha=1)
lj.pair_coeff.set('S','V',  sigma = distance_scale*2.65,epsilon = 0, alpha=1)
lj.pair_coeff.set('Au', 'Au', sigma = 1*distance_scale,epsilon = 0, alpha=1)
lj.pair_coeff.set('V', 'Au', sigma = 1*distance_scale,epsilon = 0, alpha=1)
lj.pair_coeff.set('V', 'V', sigma = 1*distance_scale,epsilon = 0, alpha=1)

##########################################
# Group particle types, set neighborlist #
##########################################

rigid = group.rigid()
nonrigid = group.nonrigid()
nlist.reset_exclusions(exclusions=['body','bond','angle','dihedral'])
nlist.set_params(check_period=1)

################
# Minimization #
################

fire= integrate.mode_minimize_rigid_fire(group=rigid, dt=0.00005)
fire_r = integrate.mode_minimize_fire(group=nonrigid,dt=0.00005)
run(5e5)
fire = integrate.mode_minimize_fire(group=nonrigid,dt=0.0005)
fire_r = integrate.mode_minimize_rigid_fire(group=rigid, dt=0.0005)
run(5e5)
fire = integrate.mode_minimize_fire(group=nonrigid,dt=0.005)
fire_r = integrate.mode_minimize_rigid_fire(group=rigid, dt=0.005)
run(5e5)
del fire
del fire_r

###################
# Zero the system #
###################

## OPTIONAL MOMENTUM ZEROER ##
#zeroer= update.zero_momentum(period=500000)
nlist.set_params(check_period=3)
integrate.mode_standard(timestep)

#########################################
# Set up logs, dump files, trajectories #
#########################################

logger = analyze.log(filename='mylog.log', period=log_period, quantities=['temperature','potential_energy','kinetic_energy','volume','pressure','pair_lj_energy','bond_harmonic_energy','angle_harmonic_energy'])

## OPTIONAL PRESSURE LOGGER ##
#pressure = analyze.log(filename='pressure.log', period=log_period, quantities=['pressure','pressure_xx','pressure_yy','pressure_zz','pressure_xy','pressure_yz','pressure_xz'])

## DUMP XML STRUCTURE FILES ##
xml = dump.xml(filename="atoms", period=xml_period)
xml.set_params(all=True)

## DUMP TRAJECTORY ##
dcd = dump.dcd(filename='mixture.dcd', period=traj_period)

######################
# Set up integrators #
######################

## OPTIONAL TEMPERATURE RAMP ##
#Tvar=variant.linear_interp([(0,t_start),(run_time,t_end)])

nonrigid_integrator=integrate.nvt(group=nonrigid, T=real_temp*kb, tau=0.65)
rigid_integrator=integrate.nvt_rigid(group=rigid, T=real_temp*kb, tau=0.65)

##############
# Box Resize #
##############
## UNCOMMENT FOR BOX RESIZING ##
#box_var=variant.linear_interp([(0,box_start),(run_time,box_end)])
#box_resize = update.box_resize(Lx=box_var,Ly=box_var,Lz=box_var)

############
# Main Run #
############
run(run_time)