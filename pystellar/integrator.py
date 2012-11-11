# -*- coding: utf-8 -*-
# 
#  integrator.py
#  pystellar
#  
#  Created by Jaberwocky on 2012-11-09.
#  Copyright 2012 Jaberwocky. All rights reserved.
# 
"""
:mod:`integrator` - RK4 Adaptive Step Size Integrator
=====================================================

"""
from __future__ import division
import numpy as np

c1,   c2,   c3,  c4,  c5,  c6 =           0,           1/4.,          3/8.,           12/13.,            1.,        1/2.
a21,                          =        1/4.,
a31, a32,                     =       3/32.,          9/32.,
a41, a42, a43,                =  1932/2197.,    -7200/2197.,    7296/2197.,
a51, a52, a53, a54            =    439/216.,            -8.,     3680/513.,      -845/4104.,
a61, a62, a63, a64, a65       =      -8/27.,             2.,   -3544/2565.,



def step_adapt(fprime,x,y,h):
    """Take an individual Runge-Kutte 4th order step."""
    x = np.asarray(x)
    y = np.asarray(y)
    assert y.ndim == 1, "Y must be a one dimensional array for a single step!"
    
    k = np.empty((6,y.shape[0]))
    
    k[0] = h * fprime(x,y)
    k[1] = h * fprime(x+a2*h,y+b21*k[0])
    k[2] = h * fprime(x+0.5*h,y+0.5*k[1])
    k[3] = h * fprime(x+h,y+k[2])
    
    ystep = y + (1/6)*(k[0] + 2*k[1] + 2*k[2] + k[3])
    return ystep
    
def step5(fprime,x,y,h):
    """docstring for step5"""
    pass