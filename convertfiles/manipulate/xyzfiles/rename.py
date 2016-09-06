import os
import sys

import manipulate_xyz as m

save='Au140.xyz'
read='Au140_sqhex.xyz'
oldname='C'
newname='Au'

m.rename_xyz(save,read,oldname,newname)

