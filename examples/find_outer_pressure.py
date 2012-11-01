#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  ps3.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-29.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

import time

import numpy as np

from astropysics.constants import Rs, Ms, Lsun
from pystellar.initial import inner_boundary, outer_boundary
from pystellar.density import mmw
from pystellar.opacity import OpacityTable
from pystellar.stellar import derivatives
from pystellar.threading import ObjectsManager,ObjectPassthrough,EngineManager,ObjectThread

X = 0.700
Y = 0.280
epsilon = 0.001
Opacity = ObjectThread(OpacityTable,ikwargs=dict(fkey='GN93hz',X=X,Y=Y),locking=True)
Opacity.start()

print "Stellar Model Photospheric Boundary Conditions"

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
start = time.clock()
ri,li,Pi,Ti = outer_boundary(R=Rs,L=Lsun,M=Ms,mu=mmw(X=X,Y=Y),optable=Opacity,Piter=10)
print "Iterable time taken: %g" % (time.clock() - start)
start = time.clock()
r,l,P,T = outer_boundary(R=Rs,L=Lsun,M=Ms,mu=mmw(X=X,Y=Y),optable=Opacity,Piter=False)
print "Absolute time taken: %g" % (time.clock() - start)
print "At m=M:"
print "    -- Algebraic -------------- | -- Iterable -------------- | -- Delta --------| "
print "    R = %10.6e cm,        | R = %10.6e cm        | dR = %11.4e |" % (r,ri,(r-ri)/r)
print "    L = %10.6e erg/s,     | L = %10.6e erg/s     | dL = %11.4e |" % (l,li,(l-li)/l)
print "    P = %10.6e Dyne/cm^2, | P = %10.6e Dyne/cm^2 | dP = %11.4e |" % (P,Pi,(P-Pi)/P)
print "    T = %10.6e K,         | T = %10.6e K         | dK = %11.4e |" % (T,Ti,(T-Ti)/T)
Opacity.stop()