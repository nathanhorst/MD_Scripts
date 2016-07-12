import nth_dihedral as nd

first='atoms.dump.0006000000.xml'
second='atoms.dump.0006050000.xml'

for i in range(11):
    x=nd.trans_v_timestep(100,str(first),str(second),i,13)
    nd.write_array(x,'trans_Tvar'+str(i)+'.txt')

