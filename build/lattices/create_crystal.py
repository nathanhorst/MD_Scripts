#!/usr/bin/env python
from numpy import *
import argparse,sys
import center as cen

examples = "EXAMPLES:\n\
create_crystal.py -s1   -a 1.00  -n 4 4 4  # simple cubic,\n\
create_crystal.py -s2   -a 1.00  -n 6 6 6  # BCC,\n\
create_crystal.py -s3   -a 1.00  -n 4 4 4  # FCC, \n\
create_crystal.py -s4   -a 1.00  -n 4 4 4  # HCP, \n\
create_crystal.py -s5   -a 1.00  -n 4 4 4 -d 'B' # AlB2, \n\
create_crystal.py -s6   -a 1.00  -n 4 4 4 -d 'B' # Li3Bi \n\ "

parser = argparse.ArgumentParser(description="create a crystal structure file",epilog=examples,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-s", dest="structure",type=int, default=1, help="choose structure (1=sc, 2=bcc, 3=fcc, 4=AlB2)")
parser.add_argument('-n', metavar='i', dest="ncells", type=int, nargs=3, default = [1, 1, 1],  help="number of unit cells in x, y and z direction (e.g. -n 3 2 4)")
parser.add_argument('-a', metavar='x', dest="lattice_constants", type=float, nargs='+', default = [1], help="supply the necessary lattice constants")
parser.add_argument('-t', dest="atomtype",type=str, default="A", help="primary atom type")
parser.add_argument('-d', dest="diatomic",type=str, help="second atom type (creates diatomic structure)")
parser.add_argument('-b', dest="basename", default="base", help="basename for output file")
parser.add_argument("-o", dest="output",type=int, default=1, help="type of output (1=xyz,2=lammps")

def main():
    args      = parser.parse_args()
    lc        = args.lattice_constants
    basename  = args.basename

    if args.structure is 1:
      print "SIMPLE CUBIC"; a=lc[0]
      cell,coords = sc(a)
    elif args.structure is 2:
      print "BCC"; a=lc[0]
      cell,coords = bcc(a)
    elif args.structure is 3:
      print "FCC"; a=lc[0]
      cell,coords = fcc(a)
    elif args.structure is 4:
      print "HCP"; a=lc[0]
      cell,coords = hcp(a)      
    elif args.structure is 5:
        print "ALB2"; a=lc[0]
        cell,coords,num = AlB2(a)
    elif args.structure is 6:
        print "Li3Bi"; a=lc[0]
        cell,coords,num = Li3Bi(a)        
    else:
      print 'structure not implemented or incorrect number of lattice constants'; return
  
    # expand crystal (n,m,k) times
    cell,coords = expand(args.ncells,cell,coords)
    

   # define atom types
    natoms = len(coords)
    types = [args.atomtype]*natoms
    if args.diatomic: 
      #for i in range(0,natoms,3): types[i] = args.diatomic
        for i in range(0,natoms):
            print(int(i)%(num[0]+num[1]))
            if int(i)%(num[0]+num[1]) >= num[0]:
                types[i] = args.diatomic
       
    # compute PBCs
    (a1,a2,a3) = cell
    a = sqrt(dot(a1,a1)); 
    b = sqrt(dot(a2,a2)); 
    c = sqrt(dot(a3,a3)); 
    if(a==0): a=1
    if(b==0): b=1
    if(c==0): c=1
    print a,b,c
    alpha = arccos(dot(a2,a3)/(b*c))*180/pi
    beta  = arccos(dot(a1,a3)/(a*c))*180/pi
    gamma = arccos(dot(a1,a2)/(a*b))*180/pi
    print dot(a1,a2),(a*b)
    print dot(a1,a2)/(a*b)
    print arccos(dot(a1,a2)/(a*b))
    print arccos(dot(a1,a2)/(a*b))*180/pi
    print '-----------------------------------------'
    print 'created ',len(coords),'atom sample'
    print '(super)cell vectors'
    print "a1 = array([% 6.2f, % 6.2f, % 6.2f])" % (a1[0],a1[1],a1[2])
    print "a2 = array([% 6.2f, % 6.2f, % 6.2f])" % (a2[0],a2[1],a2[2])
    print "a3 = array([% 6.2f, % 6.2f, % 6.2f])" % (a3[0],a3[1],a3[2])
    print '-----------------------------------------'
    print "cell for VMD : pbc set { %.3f %.3f %.3f   %.2f %.2f %.2f }; pbc box"%(a,b,c,alpha,beta,gamma)

    if   args.output is 1: writexyz(basename+".xyz",cell,coords,types)
    elif args.output is 2: writelammps(basename+".data",cell,coords)
    else: print 'output method not implemented'; return
    
    cen.center_lat_xyz()    
    
######### crystal structure functions ##########
def sc(a,r=array([0,0,0])): # simple cubic
    r0 = r
    cell = a*identity(3)
    return cell,[r0]

def bcc(a,r=array([0,0,0])): # body centered cubic
    r0 = r
    r1 = r + dot((a/2.),[1,1,1])
    cell = a*identity(3)
    return cell,[r0,r1]

def fcc(a,r=array([0,0,0])): # face centered cubic
    r0 = r
    r1 = r + dot((a/2.),[1,1,0])
    r2 = r + dot((a/2.),[1,0,1])
    r3 = r + dot((a/2.),[0,1,1])
    cell = a*identity(3)
    return cell,[r0,r1,r2,r3]

def hcp(a,r=array([0,0,0])): # hexagonal close packed
    r0 = r
    r1 = r + dot((a/2.),[1,1/sqrt(3),2*sqrt(2./3.)])
    cell = a*array([[1.0,0.0,0.0],[-0.5,0.5*sqrt(3),0.0],[0.0,0.0,1.0]])
    return cell,[r0,r1]

def AlB2(a,r=array([0,0,0])): # AlB2 base3 1A 2B
    NA=1
    NB=2    
    r0 = r
    r1 = r + dot((a/2.),[1,1/sqrt(3),1])
    r2 = r + dot((a/2.),[0,2/sqrt(3),1])
    cell = a*array([[1.0,0.0,0.0],[-0.5,0.5*sqrt(3),0.0],[0.0,0.0,1.0]])
    return cell,[r0,r1,r2],[NA,NB]

def Li3Bi(a,r=array([0,0,0])): # Li3Bi base 16 4A 12B
    NA=4
    NB=12    
    r0 = r
    r1 = r + dot((a/2.),[1,1,0])
    r2 = r + dot((a/2.),[1,0,1])
    r3 = r + dot((a/2.),[0,1,1])
    r4 = r + dot((a/2.),[0.5,0.5,0.5])
    r5 = r + dot((a/2.),[1.5,0.5,0.5])
    r6 = r + dot((a/2.),[1.5,0.5,1.5])   
    r7 = r + dot((a/2.),[0.5,0.5,1.5])
    r8 = r + dot((a/2.),[0.5,1.5,0.5])
    r9 = r + dot((a/2.),[1.5,1.5,0.5])    
    r10 = r + dot((a/2.),[1.5,1.5,1.5])
    r11 = r + dot((a/2.),[0.5,1.5,1.5])
    r12 = r + dot((a/2.),[1,0,0])    
    r13 = r + dot((a/2.),[0,1,0])
    r14 = r + dot((a/2.),[0,0,1])
    r15 = r + dot((a/2.),[1,1,1])
    cell=a*identity(3)
    return cell,[r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15],[NA,NB]
    

######## Misc. Functions #########
def expand((n1,n2,n3),cell,coords): 
    # expand current cell and coordinates by a1,a2,a3
    newcoords = []
    (a1,a2,a3) = cell
    for x in range(n1):
        for y in range(n2):
            for z in range(n3):
                r = dot(x,a1) + dot(y,a2) + dot(z,a3)
                for c in coords:
                    newcoords.append(c + r)
                    

    newcell = (dot(n1,a1),dot(n2,a2),dot(n3,a3))
    return newcell,newcoords

def periodic(cell,coords):
    # move atoms into periodic cell box
    for c in coords:
        for i in range(3):
            c[i] = c[i] % cell[i]

def translate(r,coords):
    coords = [r+c for c in coords]

# Input/Output Functions
def writexyz(file,cell,coords,types):
    # xyz
    f=open(file,'w')
    f.write(str(len(coords))+'\n')
    f.write("%12.8f %12.8f %12.8f\n" % (cell[0][0], cell[1][1], cell[2][2]))
    for i,l in enumerate(coords):
        f.write("%-2s %12.8f %12.8f %12.8f\n" % (types[i],l[0], l[1], l[2]))  # xyz format   (case for 20 unit cells, need to normalize for other values)     
        #f.write("%-2s %12.8f %12.8f %12.8f\n" % (types[i],l[0], l[1], l[2]))
        

if __name__ == "__main__":
    sys.exit(main())
