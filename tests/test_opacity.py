#!/usr/bin/env python
# 
#  test_opacity.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-04.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 


from pystellar.opacity import OpacityTable

if __name__ == '__main__':
    opacity = OpacityTable("GN93hz",load=True)
    opacity.composition(X=0.990,Y=0.000)
    opacity.setup()
    print opacity.lookup_single(logrho=0.3,logT=6.3,X=0.990,Y=0.0,Z=0.01)
    