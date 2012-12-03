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
from pystellar.density import mmw
from pystellar.opacity import OpacityTable
from pystellar.stellar import derivatives
from pystellar.threading import ObjectsManager,ObjectPassthrough,EngineManager,ObjectThread
from pystellar.star import Star

from pyshell.config import DottedConfiguration

Config = DottedConfiguration()
Config.dn = DottedConfiguration
Config.load("pystellar/Star.yml")
Config.load("Star.yml")

print "Stellar Structure Problem Set #3"
X = Config["Star.Composition.X"]
Y = Config["Star.Composition.Y"]
# P_Guess_Iteration = 2
P_Guess_Iteration = 0
Opacity = ObjectThread(OpacityTable,ikwargs=dict(fkey='GN93hz',X=X,Y=Y),locking=True)
Opacity.start()
star = Star(config=Config,optable_args=Opacity.duplicator)


print ""
print "Problem #4&5:"

print "Settings for the Sun"
print "M = %10.6g g" % Ms
print "X = %10.3f" % star.X
print "Y = %10.3f" % star.Y
print u"μ = %10.3f" % mmw(X=X,Y=Y)
print ""

print "Outer Boundary Conditions"
print "    Inital Guess:"
print "    R = %10.6e cm" % Rs
print "    L = %10.6e erg/s" % Lsun
r,l,P,T = outer_boundary(R=Rs,L=Lsun,M=Ms,mu=mmw(X=X,Y=Y),optable=Opacity,Piter=P_Guess_Iteration)
print "At m=M:"
print "    r = %10.6e cm" % r
print "    l = %10.6e erg/s" % l
print "    P = %10.6e Dyne/cm^2" % P
print "    T = %10.6e K" % T
# print "Derivatives"
# dr,dl,dP,dT = derivatives(np.array([Ms]),np.atleast_2d([r,l,P,T]),mu=mmw(X=X,Y=Y),optable=Opacity)
# print "    dr = %10.6e" % dr
# print "    dl = %10.6e" % dl
# print "    dP = %10.6e" % dP
# print "    dT = %10.6e" % dT
# print ""


# Guess from polytropes
# Pc = 1.2e17
# Tc = 1.4e7
Pc = 2.477e17
Tc = 1.571e7
m = 1e-10*Ms

print "Inner Boundary Conditions"
print "    Initial Guess: (polytrope)"
print "    P = %10.6e Dyne/cm^2" % Pc
print "    T = %10.6e K" % Tc
print "    Initial Step:"
print "    m = %10.6e g" % m
(r,l,P,T) = inner_boundary(Pc=Pc,Tc=Tc,M=Ms,mu=mmw(X=X,Y=Y),m=m,optable=Opacity,X=star.X,XCNO=star.XCNO,cfg=star.config["Data.Energy"])
print "At m=m:"
print "    r = %10.6e cm" % r
print "    l = %10.6e erg/s" % l
print "    P = %10.6e Dyne/cm^2" % P
print "    T = %10.6e K" % T
# print "Derivatives"
# (dr,dl,dP,dT) = derivatives(np.array([m]),np.atleast_2d([r,l,P,T]),mu=mmw(X=X,Y=Y),optable=Opacity,epsilon=epsilon)
# print "    dr = %10.6e" % dr
# print "    dl = %10.6e" % dl
# print "    dP = %10.6e" % dP
# print "    dT = %10.6e" % dT
# print ""
Opacity.stop()
star.stop()