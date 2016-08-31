import distance_nth_carbon as dnc
import density as d
import util as u
import radial_distribution as rd

first='atoms.dump.0006000000.xml'
second='atoms.dump.0006200000.xml'
numfiles=20
last='atoms.dump.0010800000.xml'
lattice=False
rstep=0.3
rmax=20.0
totalnp=1
nptocount=1
#length includes head and end groups
length=22

if (lattice):
    x=rd.radial_distribution(distance_array_chains(last))
    u.write_array(x,'dd.txt')

    y=dnc.dist_v_timestep(numfiles,str(first),str(second),length-1)
    u.write_array(y,'dist_end.txt')
    
    z=d.density_from_V(str(last),rstep,rmax,True)
    u.write_array(z,'density.txt')

else:
    x=rd.radial_distribution(distance_array_chains(last))
    u.write_array(x,'dd.txt')

    y=dnc.dist_v_timestep(numfiles,str(first),str(second),length-1)
    u.write_array(y,'dist_end.txt')
    
    z=d.density_from_V(str(last),rstep,rmax,True)
    u.write_array(z,'density.txt')

