from hoomd_script import *

init.create_random(N=100, phi_p=0.1)
all=group.all()
lj = pair.lj(r_cut=2.5)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
integrate.mode_standard(dt=0.005)
integrate.nvt(group=all, T=1.2, tau=0.5)
run(10000)
