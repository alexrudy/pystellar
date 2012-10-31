#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  ps3.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-29.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

from astropysics.constants import Rs, Ms, Lsun
from pystellar.initial import load_inner, load_outer
from pystellar.density import mmw
from pystellar.opacity import OpacityTable
from pystellar.threading import ObjectsManager,ObjectPassthrough

print "Stellar Structure Problem Set #2"
X = 0.700
Y = 0.280
Opacity = ObjectPassthrough(OpacityTable,nprocs=1,ikwargs=dict(fkey='GN93hz',X=X,Y=Y))
Opacity.start()

print "Problem #4:"

print "Boundaries for the Sun"
print "M = %.6g" % Ms
print "X = %.6f" % X
print "Y = %.6f" % Y
print u"μ = %.6f" % mmw(X=X,Y=Y)

print "Outer Boundary Conditions"
print "    Inital Guess:"
print "    R = %.6g" % Rs
print "    L = %.6g" % Lsun
(r,l,P,T) = load_outer(R=Rs,L=Lsun,M=Ms,mu=mmw(X=X,Y=Y),optable=Opacity)
print "At r=R:"
print "    R = %.6g" % r
print "    L = %.6g" % l
print "    P = %.6g" % P
print "    T = %.6g" % T
