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

X = 0.70
Y = 0.27

OT = ObjectsManager(OpacityTable,nprocs=1,ikwargs=dict(fkey='OP17',X=0.70,Y=0.28))
OT.start()
P = 9.987196320531581296e+04 
T = 5.564985618848681042e+03

rho = density(P=P,T=T,X=X,Y=Y)

OT.kappa(rho=rho,T=T)
print OT.retrieve()