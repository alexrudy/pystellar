#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  ps3.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-29.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 
import sys
import codecs
sys.stdout = codecs.getwriter('utf8')(sys.stdout)

import time

import numpy as np

from astropysics.constants import Rs, Ms, Lsun
from pystellar.initial import inner_boundary, outer_boundary
from pystellar.density import mmw, density
from pystellar.opacity import OpacityTable
from pystellar.stellar import derivatives
from pystellar.threading import ObjectsManager,ObjectPassthrough,EngineManager,ObjectThread
from pyshell.config import DottedConfiguration

Config = DottedConfiguration()
Config.dn = DottedConfiguration
Config.load("Star.yml")
print Config["Data.Energy"].store
X = 0.700
Y = 0.280
Opacity = ObjectThread(OpacityTable,ikwargs=dict(fkey='GN93hz',X=X,Y=Y,snap=True),locking=True)
Opacity.start()

print ""

print "Settings for the Sun"
print "M = %10.6g g" % Ms
print "X = %10.3f" % X
print "Y = %10.3f" % Y
print u"μ = %10.3f" % mmw(X=X,Y=Y)
print ""

m = 1.9890999983295044e+33

print "Inner Boundary Conditions"
print "    Initial Step:"
print "    m = %10.6e g" % m
(r,l,P,T) = [  6.96000000e+10,   2.13492744e+33,   3.64246107e+27,  5.77489402e+03]
print "At m=m:"
print "    r = %10.6e cm" % r
print "    l = %10.6e erg/s" % l
print "    P = %10.6e Dyne/cm^2" % P
print "    T = %10.6e K" % T
print u"      = %10.6e " % density(P=P,T=T,X=X,Y=Y)
print "Derivatives"
(dr,dl,dP,dT) = derivatives(np.array([m]),np.atleast_2d([r,l,P,T]),mu=mmw(X=X,Y=Y),optable=Opacity,X=X,XCNO=(1-X-Y),cfg=Config["Data.Energy"])
print "    dr = %10.6e" % dr
print "    dl = %10.6e" % dl
print "    dP = %10.6e" % dP
print "    dT = %10.6e" % dT
print ""
Opacity.stop()