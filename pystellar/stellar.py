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

import numpy as np

from .density import density, ldensity
from .energy import epp, eCNO

def drdm(r,rho):
    ur"""Find the derivative of radius with respect to mass.
    
    .. math::
        \frac{dr}{dm} = \frac{1}{4\pi r^2 \rho}
        
    
    :param r: Radius
    :param rho: Density
    
    """
    return 1 / (4*np.pi*np.power(r,2)*rho)
    
def dlogrdlogm(lr,lrho):
    ur"""Find the derivative of radius with respect to mass.
    
    .. math::
        \frac{d\log(r)}{d\log(m)} = \frac{1}{4\pi r^2 \rho}
        
    
    :param r: Radius
    :param rho: Density
    
    """
    return np.log10(1 / (4*np.pi)) - 2 * lr - lrho
    
def dPdm(r,m):
    ur"""Find the derivative of the pressure with respect to mass.
    
    .. math::
        \frac{dP}{dm} = -\frac{Gm}{4\pi r^2}
        
    :param r: Radius enclosed
    :param m: Mass enclosed
    """
    from .constants import G
    return (- G * m) / (4*np.pi*np.power(r,4))
    
def dldm(T,rho,X,XCNO,cfg):
    ur"""Find the derivative of luminosity with respect to mass.
    
    """
    return epp(T=T,rho=rho,X=X,c=cfg) + eCNO(T=T,rho=rho,X=X,XCNO=XCNO,c=cfg)

def dlogldlogm(T,rho,X,XCNO,cfg):
    """Logarithmic Luminosity Derivative"""
    return np.log10(dldm(T,rho,X,XCNO,cfg))

def dTdm(m,r,l,P,T,rho,optable):
    ur"""Derivative of tempertature with respect to mass.
    
    .. math::
        
        \frac{dT}{dm} = - \frac{GmT}{4\pi r^4 P} \nabla
        
    :param m: Mass enclosed
    :param T: Temperature
    :param r: Radius enclosed
    :param P: Pressure enclosed
    :param grad: The temperature gradient locally.
    
    .. todo::
        Figure out how this function should actually handle the temperature gradient
    """
    
    from .constants import G, gradT_ad
    rgrad = radiative_gradient(T=T,P=P,l=l,m=m,rho=rho,optable=optable)
    grad = np.zeros(rgrad.shape)
    grad[:] = rgrad # Schwartzchild criterion
    grad[rgrad < gradT_ad] = gradT_ad
    return -(G*m*T)/(4 * np.pi * np.power(r,4) * P) * grad
    
def radiative_gradient(T,P,l,m,rho,optable):
    """docstring for radiative_gradient"""
    from .constants import G, c, a
    k = optable.retrieve()    
    return 3/(16*np.pi*a*c*G) * (k * l * P)/(m * np.power(T,4))
    
    
def derivatives(xs,ys,mu,optable,X,XCNO,cfg):
    """Return the full set of derivatives.
    
    :param xs: The mass positions, m, at which to find the derivative.
    :param ys: The sets of (r,l,P,T) at which to find the derivative.
    :returns: The sets of (dr,dl,dP,dT) of derivatives.
    """
    
    r, l, P, T = np.asarray(ys).T        
    m = np.asarray(xs)
    rho = density(P=P,T=T,mu=mu)
    optable.kappa(rho=rho,T=T)
    dr = drdm(r=r,rho=rho)
    dl = dldm(T=T,rho=rho,X=X,XCNO=XCNO,cfg=cfg)
    dP = dPdm(r=r,m=m)
    dT = dTdm(m=m,T=T,r=r,P=P,l=l,rho=rho,optable=optable)
    return np.vstack((dr,dl,dP,dT))
    
def log_derivatives(xs,ys,mu,optable,X,XCNO,cfg):
    """Return the full set of derivatives.
    
    :param xs: The mass positions, m, at which to find the derivative.
    :param ys: The sets of (\log(r),\log(l),P,T) at which to find the derivative.
    :returns: The sets of (dr,dl,dP,dT) of derivatives.
    """
    r, l, P, T = np.asarray(ys).T
    lm = np.asarray(xs)
    m = np.power(10,lm)
    logspace = m * np.log(10)
    rho = density(P=P,T=T,mu=mu)
    lrho = ldensity(logP=np.log10(P),logT=np.log10(T),mu=mu)
    optable.kappa(rho=rho,T=T)
    dr = drdm(r=r,rho=rho)
    dl = dldm(T=T,rho=rho,X=X,XCNO=XCNO,cfg=cfg)
    dP = dPdm(r=r,m=m)
    dT = dTdm(m=m,T=T,r=r,P=P,l=l,rho=rho,optable=optable)
    return np.vstack((dr,dl,dP,dT)) * logspace
    
    