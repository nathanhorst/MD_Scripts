from __future__ import division
import math
import sys
import os
import manipulate_xml as man

save='trans_shift.xml'
read='trans_chain.xml'
x=-18
y=0
z=0
npart=32


man.shift_xml(save,read,x,y,z,npart)
