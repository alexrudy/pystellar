#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  ps3.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-29.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

import numpy as np

from astropysics.constants import Rs, Ms, Lsun
from pystellar.initial import inner_boundary, outer_boundary
from pystellar.density import mmw
from pystellar.opacity import OpacityTable
from pystellar.stellar import derivatives
from pystellar.threading import ObjectsManager,ObjectPassthrough,EngineManager

print "Stellar Structure Problem Set #2"
X = 0.700
Y = 0.280
epsilon = 0.001
Opacity = ObjectsManager(OpacityTable,nprocs=1,ikwargs=dict(fkey='GN93hz',X=X,Y=Y))
Opacity.start()

print "Problem #4:"

print "Settings for the Sun"
print "M = %10.6g g" % Ms
print "X = %10.3f" % X
print "Y = %10.3f" % Y
print u"ε = %10.3f ergs/g" % epsilon
print u"μ = %10.3f" % mmw(X=X,Y=Y)

print "Outer Boundary Conditions"
print "    Inital Guess:"
print "    R = %10.6e cm" % Rs
print "    L = %10.6e erg/s" % Lsun
(r,l,P,T) = outer_boundary(R=Rs,L=Lsun,M=Ms,mu=mmw(X=X,Y=Y),optable=Opacity)
print "At m=M:"
print "    R = %10.6e cm" % r
print "    L = %10.6e erg/s" % l
print "    P = %10.6e Dyne/cm^2" % P
print "    T = %10.6e K" % T
print "Derivatives"
(dr,dl,dP,dT) = derivatives(np.array([Ms]),np.atleast_2d([r,l,P,T]),mu=mmw(X=X,Y=Y),optable=Opacity,epsilon=epsilon)
print "    dr = %10.6e" % dr
print "    dl = %10.6e" % dl
print "    dP = %10.6e" % dP
print "    dT = %10.6e" % dT
print ""


# Guess from polytropes
Pc = 1.2e17
Tc = 1.4e7
m = 1e-30*Ms

print "Inner Boundary Conditions"
print "    Initial Guess: (polytrope)"
print "    P = %10.6e Dyne/cm^2" % Pc
print "    T = %10.6e K" % Tc
print "    Initial Step:"
print "    m = %10.6e g" % m
(r,l,P,T) = inner_boundary(Pc=Pc,Tc=Tc,M=Ms,mu=mmw(X=X,Y=Y),m=m,epsilon=epsilon,optable=Opacity)
print "At M=m:"
print "    R = %10.6e cm" % r
print "    L = %10.6e erg/s" % l
print "    P = %10.6e Dyne/cm^2" % P
print "    T = %10.6e K" % T
print "Derivatives"
(dr,dl,dP,dT) = derivatives(np.array([m]),np.atleast_2d([r,l,P,T]),mu=mmw(X=X,Y=Y),optable=Opacity,epsilon=epsilon)
print "    dr = %10.6e" % dr
print "    dl = %10.6e" % dl
print "    dP = %10.6e" % dP
print "    dT = %10.6e" % dT
print ""
Opacity.stop()