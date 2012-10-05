#!/usr/bin/env python
# 
#  read_opacity.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-03.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

from pystellar.opacity import read_table_opal, get_interpolator

if __name__ == '__main__':
    interp = get_interpolator(*read_table_opal("GN93hz"))
    print interp((2,5,1.0,0.0,0.0,0.0,0.0))
    