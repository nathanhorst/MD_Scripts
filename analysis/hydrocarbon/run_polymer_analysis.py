import nth_dihedral as nd
import distance_nth_carbon as dnc
import density as d
import util as u

first='atoms.dump.0000300000.xml'
second='atoms.dump.0000500000.xml'
numfiles=5
last='atoms.dump.0001100000.xml'
lattice=True
rstep=0.8/3.96
rmax=6.0
totalnp=32
nptocount=1
#length includes head and end groups
length=13

if (lattice):
    for i in range(length-2):
        x=nd.trans_v_timestep_lat(numfiles,str(first),str(second),i,length,totalnp,nptocount)
        u.write_array(x,'trans_Tvar'+str(i)+'.txt')
    #print i


### dist_v_timestep(#dumpfiles,firstfile,secondfile,Lpoly-1) ###

    y=dnc.dist_v_timestep(numfiles,str(first),str(second),length-1)
    u.write_array(y,'dist_12thC.txt')
    
    z=d.density_from_V(str(last),rstep,rmax,True)
    u.write_array(z,'density.txt')

#a=d.histogram(d.inter_angle(last,12),30)
#u.write_array(a,'angle12thC.txt')

else:
    for i in range(length-2):
        x=nd.trans_v_timestep(numfiles,str(first),str(second),i,length)
        u.write_array(x,'trans_Tvar'+str(i)+'.txt')


### dist_v_timestep(#dumpfiles,firstfile,secondfile,Lpoly-1) ###

    y=dnc.dist_v_timestep(numfiles,str(first),str(second),length-1)
    u.write_array(y,'dist_12thC.txt')
    
    z=d.density_from_V(str(last),rstep,rmax,True)
    u.write_array(z,'density.txt')

#a=d.histogram(d.inter_angle(last,12),30)
#u.write_array(a,'angle12thC.txt')
