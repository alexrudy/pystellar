#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  opacity.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-10.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

from pystellar.opacity import OpacityTable

print "\nProblem #6:"

logTs   = [6.3, 5.0]
logrhos = [0.3,-4.0]
X = 0.700
Y = 0.280
opacity = OpacityTable("GN93hz",load=False)
opacity.composition(X=X,Y=Y)
print u"Fixed composition at X=%.3g,Y=%.3g,Z=%.3g" % (opacity.X,opacity.Y,opacity.Z)
print u"Using table %d" % opacity.n

for logT,logrho in zip(logTs,logrhos):
    kappa = opacity.lookup(logrho=logrho,logT=logT)
    print u"log(T)=%5.3f, log(ρ)=%6.3f\n → κ=%5.3f" % (logT,logrho,kappa)
    
opacity.composition(X=0.00,Y=1.00)
logrho, logT = opacity.invert_points(logR=-8.0,logT=3.75).T
kappa = opacity.lookup(logrho=logrho,logT=logT)
print u"log(T)=%5.3f, log(ρ)=%6.3f\n → κ=%5.3f" % (logT,logrho,kappa)