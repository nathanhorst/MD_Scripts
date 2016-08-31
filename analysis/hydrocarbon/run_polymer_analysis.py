import nth_dihedral as nd
import distance_nth_carbon as dnc
import density as d
import util as u
import intermolecular_angle as ia

lattice=False

numfiles=50
first='atoms.dump.0004500000.xml'
second='atoms.dump.0004600000.xml'
totalnp=32
nptocount=1
#length includes head and end groups
length=13

last='atoms.dump.0009400000.xml'
rstep=0.5/3.96
rmax=6.0
volume_division=True

boxes=50

if (lattice):
    
    x=nd.trans_v_timestep_lat_matrix(numfiles,first,second, length,totalnp,nptocount)
    y=nd.n_files(numfiles,first,second)
    for i in range(0,length-2):
        u.write_array(nd.trans_v_timestep_combine(y,i,x,length),'trans_Tvar'+str(i)+'.txt') 
    
    x= dnc.dist_v_timestep_lat_matrix(numfiles,first,second,length,totalnp,nptocount)
    y= u.n_files(numfiles,first,second)    
    for i in range(1,length):
        u.write_array(dnc.nth_distance_lat(y,i,x,length),'dist'+str(i)+'thC.txt')
    
    z=d.density_from_V(str(last),rstep,rmax,volume_division)
    u.write_array(z,'density.txt')

   

else:
    x=nd.trans_v_timestep_matrix(numfiles,first,second, length)
    y=nd.n_files(numfiles,first,second)
    for i in range(0,length-2):
        u.write_array(nd.trans_v_timestep_combine(y,i,x,length),'trans_Tvar'+str(i)+'.txt') 

    x= dnc.dist_v_timestep_lat_matrix(numfiles,first,second,length,1,1)
    y= nd.n_files(numfiles,first,second)    
    for i in range(1,length):
        u.write_array(dnc.nth_distance_lat(y,i,x,length),'dist'+str(i)+'thC.txt')
    
    
    z=d.density_from_V(str(last),rstep,rmax,volume_division)
    u.write_array(z,'density.txt')

    a=u.histogram(ia.inter_angle(last,length-1),boxes)
    u.write_array(a,'angle12thC.txt')
