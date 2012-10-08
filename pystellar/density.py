# 
#  density.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-07.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 
import numpy as np
import scipy as sp

from pkg_resources import resource_filename
from warnings import warn

import logging

log = logging.getLogger(__name__)

def mmw(X,Y):
    """docstring for mean_mol_weight"""
    mmw = (2.0) / (1.0 + 3.0 * X + 0.5 * Y)
    return mmw
    
def ldensity(logP,logT,mu):
    """docstring for ldensity"""
    P = 10.0**logP
    T = 10.0**logT
    return density(P,T,mu)
    
def density(P,T,mu):
    """Density for a given Pressure, Temperature, and mean mol. weight"""
    from .constants import a, mh, kb
    rho = (mu * mh) / (kb * T) * (P - (a/3)*T**4)
    return rho
    
def lbeta(logP,logT,mu):
    """docstring for lbeta"""
    P = 10.0**logP
    T = 10.0**logT
    return beta(P,T,mu)
    
def beta(P,T,mu):
    """Ratio of gas to total pressure"""
    from .constants import a, mh, kb
    Prad = (a/3) * T**4
    beta = (1 - (Prad/P))
    return beta