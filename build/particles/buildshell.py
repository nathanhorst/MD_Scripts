# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 09:16:36 2016

@author: nathan
"""

from __future__ import division
 
import numpy as np

import nanoparticle_core as npc


uc_size=2.5
uc_shellsize=uc_size-0.5
uc_bcc=npc.NanoBcc(uc_size)
uc_shell=npc.NanoBcc(uc_shellsize)
#uc_fcc = npc.NanoFcc(uc_size)

searchquery1='C'
vec_all=[]
TYPE='Au'

fin= open("Au201shell.xyz",'w')

with open('base.xyz') as f1:
    lines=f1.readlines()
    for i, line in enumerate(lines):
        if line.startswith(searchquery1):
            split=str(line).split()
            vec_all.append(np.array([split[1],split[2],split[3]]).astype(np.float))
    V=np.array([[0.0,0.0,0.0]])
    count=1
    for ind, vec in enumerate(vec_all):
        if uc_bcc.check_point(np.array(vec)):
            if not uc_shell.check_point(np.array(vec)):
                count+=1
                V=np.append(V,[vec],axis=0)

print count
fin.write(('%s \n\n')%(count))
for i in range(1,count):
    fin.write(('%s %1.6f %1.6f %1.6f \n')%(TYPE, V[i][0], V[i][1], V[i][2]))
fin.write(('V %1.6f %1.6f %1.6f \n')%(V[0][0], V[0][1], V[0][2]))
           

