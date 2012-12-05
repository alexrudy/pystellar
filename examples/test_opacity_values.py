#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  test_opacity_values.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-12-02.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 


from __future__ import division

from pystellar.density import density
from pystellar.opacity import OpacityTable
from pystellar.threading import ObjectsManager, EngineManager

import numpy as np
import time
import logging

X = 0.70
Y = 0.28

log = logging.getLogger('pystellar.opacity')
log.setLevel(logging.DEBUG)
console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
log.addHandler(console)

OP = OpacityTable(fkey='OP17',X=0.70,Y=0.28,efkey='cunha06')
P, T = [8.726086186677213013e+07, 4.576702504411734481e+03]

rho = density(P=P,T=T,X=X,Y=Y)

print OP.kappa(rho=rho,T=T)
print OP.make_points(logT=np.log10(T),logrho=np.log10(rho))
T= 11264.2381423
rho=0.0044603404639
print OP.kappa(T=T,rho=rho)