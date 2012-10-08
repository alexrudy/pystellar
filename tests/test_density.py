#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  test_density.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-07.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

from __future__ import division

from pystellar.density import ldensity, mmw

from pystellar.constants import a, mh, kb

print "a=%5.3g" % a
print "Prad = %5.3g" % (a/3.0 * (10**(6.91))**4.0)
print "Pgas = %5.3g" % (10**16.87 - (a/3.0 * (10**(6.91))**4.0))
print "mu = %5.3g" % mmw(X=0.7,Y=0.28)