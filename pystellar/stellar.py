# -*- coding: utf-8 -*-
# 
#  stellar.py
#  pystellar
#  Equations of stellar state.
#  
#  Created by Alexander Rudy on 2012-10-29.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 
"""
:mod:`stellar` - Equations of State
===================================

"""
from __future__ import division

from .density import density

def drdm(r,rho):
    ur"""Find the derivative of radius with respect to mass.
    
    .. math::
        \frac{dr}{dm} = \frac{1}{4\pi r^2 \rho}
        
    
    :param r: Radius
    :param rho: Density
    
    """
    return 1 / (4*np.pi*np.power(r,2),rho)
    
def dPdm(r,m):
    ur"""Find the derivative of the pressure with respect to mass.
    
    .. math::
        \frac{dP}{dm} = -\frac{Gm}{4\pi r^2}
        
    :param r: Radius enclosed
    :param m: Mass enclosed
    """
    from .constants import G
    return (- G * m) / (4*np.pi*np.power(r,2))
    
def dldm(epsilon):
    ur"""Find the derivative of luminosity with respect to mass.
    
    .. todo::
        This function is entirely dependent on the total energy generation rate, a quantity we haven't explored yet at all. As such, this function really does nothing!
    """
    return epsilon

def dTdm(m,T,r,P,grad):
    ur"""Derivative of tempertature with respect to mass.
    
    .. math::
        
        \frac{dT}{dm} = - \frac{GmT}{4\pi r^4 P} \nabla
        
    :param m: Mass enclosed
    :param T: Temperature
    :param r: Radius enclosed
    :param P: Pressure enclosed
    :param grad: The temperature gradient locally.
    """
    from .constants import G
    return -(G*m*T)/(4 * np.pi * np.power(r,4) * P)
    
def derivatives(xs,ys,mu,epsilon,grad):
    """Return the full set of derivatives.
    
    :param xs: The mass positions, m, at which to find the derivative.
    :param ys: The sets of (r,l,P,T) at which to find the derivative.
    :returns: The sets of (dr,dl,dP,dT) of derivatives.
    """
    
    r, l, P, T = np.asarray(ys).T
    m = np.asarray(xs)
    rho = density(P=P,T=T,mu=mu)
    dr = drdm(r=r,rho=rho)
    dl = dldm(epsilon)
    dP = dPdm(r=r,m=m)
    dT = dTdm(m=m,T=T,r=r,P=P,grad=grad)
    return np.hstack((dr,dl,dP,dT)).T
    