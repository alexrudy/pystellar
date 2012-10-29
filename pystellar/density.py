# 
#  density.py
#  Calculation of Density from logT, logP, X & Y
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-07.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 
"""
:mod:`density` - Density Calculation
====================================

Density & Beta
**************

.. autofunction::
    density
    
.. autofunction::
    beta

Logarithmic Functions
*********************

.. autofunction::
    ldensity
    
.. autofunction::
    lbeta

Support Functions
*****************

.. autofunction::
    mmw

"""
from __future__ import division
import numpy as np
import scipy as sp

from pkg_resources import resource_filename
from warnings import warn

import logging

log = logging.getLogger(__name__)

def mmw(X,Y):
    """Calculate the Mean Molecular Weight"""
    mmw = (2.0) / (1.0 + 3.0 * X + 0.5 * Y)
    return mmw
    
def ldensity(logP,logT,mu=None,X=None,Y=None):
    """Density for a given log(P), log(T) and mean molecular weight."""
    P = np.power(10,logP)
    T = np.power(10,logT)
    return density(P,T,mu,X,Y)
    
def density(P,T,mu=None,X=None,Y=None):
    """Density for a given Pressure, Temperature, and mean mol. weight."""
    from .constants import a, mh, kb
    if mu is None:
        mu = mmw(X,Y)
    rho = (mu * mh) / (kb * T) * (P - (a/3) * np.power(T,4))
    return rho
    
def lbeta(logP,logT,mu=None,X=None,Y=None):
    """Calculate the pressure ration (beta) from log(P), log(T) and mean molecular weight."""
    P = np.power(10,logP)
    T = np.power(10,logT)
    return beta(P,T,mu,X,Y)
    
def beta(P,T,mu=None,X=None,Y=None):
    """Ratio of gas to total pressure"""
    from .constants import a, mh, kb
    if mu is None:
        mu = mmw(X,Y)
    Prad = (a/3) * np.power(T,4)
    beta = (1 - (Prad/P))
    return beta