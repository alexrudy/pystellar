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

print "1: OP17"
print "2: GN93hz"
logTs   = [6.3, 5.0]
logrhos = [0.3,-4.0]
X = 0.700
Y = 0.280
opacity1 = OpacityTable("OP17",load=False)
opacity2 = OpacityTable("GN93hz",load=False)
opacity1.composition(X=X,Y=Y)
opacity2.composition(X=X,Y=Y)
print u"1: Fixed composition at X=%.3g,Y=%.3g,Z=%.3g" % (opacity1.X,opacity1.Y,opacity1.Z)
print u"2: Fixed composition at X=%.3g,Y=%.3g,Z=%.3g" % (opacity2.X,opacity2.Y,opacity2.Z)
print u"1: Using table %d" % opacity1.n
print u"2: Using table %d" % opacity2.n

for logT,logrho in zip(logTs,logrhos):
    kappa = opacity1.lookup(logrho=logrho,logT=logT)
    print u"1: log(T)=%5.3f, log(ρ)=%6.3f → κ=%5.3f" % (logT,logrho,kappa)
    kappa = opacity2.lookup(logrho=logrho,logT=logT)
    print u"2: log(T)=%5.3f, log(ρ)=%6.3f → κ=%5.3f" % (logT,logrho,kappa)
    
    
opacity1.composition(X=0.00,Y=1.00)
opacity2.composition(X=0.00,Y=1.00)
print u"1: Fixed composition at X=%.3g,Y=%.3g,Z=%.3g" % (opacity1.X,opacity1.Y,opacity1.Z)
print u"2: Fixed composition at X=%.3g,Y=%.3g,Z=%.3g" % (opacity2.X,opacity2.Y,opacity2.Z)
print u"1: Using table %d" % opacity1.n
print u"2: Using table %d" % opacity2.n

logrho, logT = opacity1.invert_points(logR=-8.0,logT=3.75).T
try:
    kappa = opacity1.lookup(logrho=logrho,logT=logT)
    print u"1: log(T)=%5.3f, log(ρ)=%6.3f → κ=%5.3f" % (logT,logrho,kappa)
except ValueError:
    print u"1: log(T)=%5.3f, log(ρ)=%6.3f → κ=NaN" % (logT,logrho)
    
try:
    kappa = opacity2.lookup(logrho=logrho,logT=logT)
    print u"2: log(T)=%5.3f, log(ρ)=%6.3f → κ=%5.3f" % (logT,logrho,kappa)
except ValueError:
    print u"2: log(T)=%5.3f, log(ρ)=%6.3f → κ=NaN" % (logT,logrho)