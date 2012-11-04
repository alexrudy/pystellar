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
print "Testing Outer Presure Conversion Mehtods"
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
print u"    R = %10.6e cm,        | R = %10.6e cm        | ΔR = %11.4e |" % (r,ri,(r-ri)/r)
print u"    L = %10.6e erg/s,     | L = %10.6e erg/s     | ΔL = %11.4e |" % (l,li,(l-li)/l)
print u"    P = %10.6e Dyne/cm^2, | P = %10.6e Dyne/cm^2 | ΔP = %11.4e |" % (P,Pi,(P-Pi)/P)
print u"    T = %10.6e K,         | T = %10.6e K         | ΔT = %11.4e |" % (T,Ti,(T-Ti)/T)


Pc = 2.477e17
Tc = 1.571e7
m = 1e-30*Ms
print ""
print "Inner Boundary Conditions"
print "Testing Delta Movement of Temperature & Pressure from Boundary"
print "    Initial Guess: (polytrope)"
print "    P = %10.6e Dyne/cm^2" % Pc
print "    T = %10.6e K" % Tc
print "    Initial Step:"
print "    m = %10.6e g" % m
(rs,ls,Ps,Ts) = inner_boundary(Pc=Pc,Tc=Tc,M=Ms,mu=mmw(X=X,Y=Y),m=m,epsilon=epsilon,optable=Opacity)
(rc,lc,Pc,Tc) = (0,0,Pc,Tc)
print "At m=m:"
print "    r = %10.6e cm" % rs
print "    l = %10.6e erg/s" % ls
print "    P = %10.6e Dyne/cm^2" % Ps
print "    T = %10.6e K" % Ts
print "    -- Initial ---------------- | -- 1st Step -------------- | -- Delta --------| "
print u"    R = %10.6e cm,        | R = %10.6e cm        | ΔR = %11.4e |" % (rs,rc,abs(rs-rc)/rs)
print u"    L = %10.6e erg/s,     | L = %10.6e erg/s     | ΔL = %11.4e |" % (ls,lc,abs(ls-lc)/ls)
print u"    P = %10.6e Dyne/cm^2, | P = %10.6e Dyne/cm^2 | ΔP = %11.4e |" % (Ps,Pc,abs(Ps-Pc)/Ps)
print u"    T = %10.6e K,         | T = %10.6e K         | ΔT = %11.4e |" % (Ts,Tc,abs(Ts-Tc)/Ts)




print "Derivatives"
print "Testing in Array Mode"
dT, dP, dl, dr = derivatives(np.array([m]),np.atleast_2d([[ri,li,Pi,Ti],[rs,ls,Ps,Ts]]),mu=mmw(X=X,Y=Y),optable=Opacity,epsilon=epsilon)
# print derivs
print "    dr = [%14.6e,%14.6e]" % tuple(dr)
print "    dl = [%14.6e,%14.6e]" % tuple(dl)
print "    dP = [%14.6e,%14.6e]" % tuple(dP)
print "    dT = [%14.6e,%14.6e]" % tuple(dT)
# print ""



Opacity.stop()