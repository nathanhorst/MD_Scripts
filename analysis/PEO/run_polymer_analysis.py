#import nth_dihedral as nd
import distance_nth_carbon as dnc
import density as d
import util as u

first='atoms.dump.0000500000.xml'
second='atoms.dump.0000552083.xml'

last='atoms.dump.0010499936.xml'

rstep=0.8
rmax=20.0

y=dnc.dist_v_timestep(100,str(first),str(second),22)
u.write_array(y,'dist_12thC.txt')

z=d.density(str(last),rstep,rmax,True)
u.write_array(z,'density.txt')

#a=d.histogram(d.inter_angle(last,12),30)
#u.write_array(a,'angle12thC.txt')

