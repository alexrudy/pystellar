#!/usr/bin/env python
# 
#  opacity.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-04.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 


from pystellar.opacity import OpacityTable

if __name__ == '__main__':
    opacity = OpacityTable("GN93hz")
    print opacity.lookup_single(logrho=0.3,logT=6.3,X=0.1,Y=0.9,Z=0.0)
    