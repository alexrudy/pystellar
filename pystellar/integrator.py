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
a51, a52, a53, a54,           =    439/216.,            -8.,     3680/513.,      -845/4104.,
a61, a62, a63, a64, a65,      =      -8/27.,             2.,   -3544/2565.,      1859/4104.,          -11/40,
b1 ,  b2,  b3,  b4,  b5,   b6 =     16/135.,             0.,   6656/12825.,    28561/56430.,           -9/50,      2/55.
b1s, b2s, b3s, b4s, b5s,  b6s =     25/216.,             0.,    1408/2565.,      2197/4104.,            -1/5,         0


def step_adapt(f,x,y,h,args=(),tol=1e-8):
    """Take an individual Runge-Kutte 4th order step."""
    x = np.asarray(x)
    y = np.asarray(y)
    assert y.ndim == 1, "Y must be a one dimensional array for a single step!"
    fprime = lambda x,y : f(y,x,*args)
    k = np.empty((6,y.shape[0]))
    
    k[0] = h * fprime(x+c1,y)
    k[1] = h * fprime(x+c2*h,y+a21*k[0])
    k[2] = h * fprime(x+c3*h,y+a31*k[0]+a32*k[1])
    k[3] = h * fprime(x+c4*h,y+a41*k[0]+a42*k[1]+a43*k[2])
    k[4] = h * fprime(x+c5*h,y+a51*k[0]+a52*k[1]+a53*k[2]+a54*k[3])
    k[5] = h * fprime(x+c6*h,y+a61*k[0]+a62*k[1]+a63*k[2]+a64*k[3]+a65*k[4])
    
    ystep = y + (b1s*k[0] + b2s*k[1] + b3s*k[2] + b4s*k[3] + b5s*k[4] + b6s*k[5])
    err = (b1*k[0] + b2*k[1] + b3*k[2] + b4*k[3] + b5*k[4] + b6*k[5]) - (b1s*k[0] + b2s*k[1] + b3s*k[2] + b4s*k[3] + b5s*k[4] + b6s*k[5])
    hn = 0.9 * h * tol/np.abs(err)
    return ystep,h
    
    
def integrate(fprime,x,y0,h0,tol=1e-8):
    """Integrator!"""
    x = np.asarray(x)
    y0 = np.asarray(y0)
    y = np.array((x.shape[0],y0.shape[0]))
    hc = h0
    xc = x[0]
    xn = xc
    yc = y0
    while xc <= x[-1]:
        xn += hc
        yc,hc = step_adapt(fprime,xc,yc,hc,tol=tol)
        xc = xn
        # y[x == xc] = yc
    return x,y,xc,yc
    