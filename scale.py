import os
import sys

import manipulate_xyz as m

save='complete_down.xyz'
read='complete.xyz'
shiftx=1.0/3.96
shifty=1.0/3.96
shiftz=1.0/3.96

m.scale_xyz(save,read,shiftx,shifty,shiftz)

