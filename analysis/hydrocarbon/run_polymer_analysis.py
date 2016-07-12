import nth_dihedral as nd
import distance_nth_carbon as dnc
import density as d

first='atoms.dump.0006000000.xml'
second='atoms.dump.0006050000.xml'

last='atoms.dump.0010950000.xml'

rstep=0.8/3.96
rmax=6.0

for i in range(11):
    x=nd.trans_v_timestep(100,str(first),str(second),i,13)
    nd.write_array(x,'trans_Tvar'+str(i)+'.txt')

y=dnc.dist_v_timestep(100,str(first),str(second),12)
nd.write_array(y,'dist_12thC.txt')

z=d.density(str(last),rstep,rmax)
nd.write_array(z,'density.txt')

#a=d.histogram(d.inter_angle(last,12),30)
#nd.write_array(a,'angle12thC.txt')

