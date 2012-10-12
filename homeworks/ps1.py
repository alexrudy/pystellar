#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  ps1.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-07.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

from pystellar.density import ldensity, lbeta, mmw
from pystellar.opacity import OpacityTable
import numpy as np

print "Stellar Structure Problem Set 1"

# Problem #5
print "\nProblem #5:"

Xs    = [0,     0.70]
Ys    = [0.98,  0.28]
logTs = [7.55,  6.91]
logPs = [16.85,16.87]


for X,Y,logT,logP in zip(Xs,Ys,logTs,logPs):
    mu = mmw(X=X,Y=Y)
    beta = lbeta(logT=logT,logP=logP,mu=mu)
    dens = ldensity(logT=logT,logP=logP,mu=mu)
    print u"For X=%5.3f, Y=%5.3f, log(T)=%5.3f, log(P)=%5.3f\n → μ=%5.3f → ρ=%5.3f, β=%5.3f" % (X,Y,logT,logP,mu,dens,beta)


# Problem #6
print "\nProblem #6:"

logTs   = [6.3, 5.0]
logrhos = [0.3,-4.0]
X = 0.700
Y = 0.280
opacity = OpacityTable("GN93hz",load=True)
opacity.composition(X=X,Y=Y)
print u"Fixed composition at X=%.3g,Y=%.3g,Z=%.3g" % (opacity.X,opacity.Y,opacity.Z)
print u"Using table %d" % (opacity.n + 1)

for logT,logrho in zip(logTs,logrhos):
    kappa = opacity.lookup(logrho=logrho,logT=logT)
    print u"log(T)=%5.3f, log(ρ)=%6.3f\n → log(κ)=%5.3f" % (logT,logrho,kappa)
