# -*- coding: utf-8 -*-
# 
#  initial.py
#  pystellar
#  
#  Created by Jaberwocky on 2012-10-29.
#  Copyright 2012 Jaberwocky. All rights reserved.
# 
"""
:mod:`initial` - Initial conditions
===================================

Boundary Condition Functions
****************************

.. autofunction::
    load_inner

.. autofunction::
    load_outer

Support Functions
*****************
    
Inner Boundary
~~~~~~~~~~~~~~
    
.. autofunction::
    inner_radius

.. autofunction::
    inner_pressure
    
.. autofunction::
    inner_temperature

Outer Boundary
~~~~~~~~~~~~~~

.. autofunction::
    outer_temperature
    
.. autofunction::
    outer_pressure

"""


from __future__ import division
import numpy as np

from .density import density, ldensity
from .stellar import radiative_gradient
from .energy import epp, eCNO

import logging

log = logging.getLogger(__name__)

def inner_boundary(Pc,Tc,M,mu,m,optable,X,XCNO,cfg,convective=False):
    r"""Load our inner boundary conditions for our star.
    
    :param Pc: Central Pressure, :math:`P_c`
    :param Tc: Central Temperature, :math:`T_c`
    :param M: Total Star Mass
    :param mu: Mean Molecular Weight, :math:`\mu`
    :param m: Initial mass step from center of star.
    :returns: Tuple of `(r,l,P,T)`, the four equation of state values at point m.
        
    """
    rhoc = density(P=Pc,T=Tc,mu=mu)
    epsilon = (epp(T=Tc,rho=rhoc,X=X,c=cfg) + eCNO(T=Tc,rho=rhoc,X=X,XCNO=XCNO,c=cfg))
    r = inner_radius(rho=rhoc,m=m)
    l = m * epsilon
    P = inner_pressure(Pc=Pc,rho=rhoc,m=m)
    T = inner_temperature(Tc=Tc,Pc=Pc,rho=rhoc,m=m,epsilon=epsilon,optable=optable,convective=convective)
    return (r,l,P,T)
    
def log_inner_boundary(Pc,Tc,lM,mu,lm,optable,X,XCNO,cfg,convective=True):
    r"""Load our inner boundary conditions for our star in Log Space.
    
    :param Pc: Central Pressure, :math:`P_c`
    :param Tc: Central Temperature, :math:`T_c`
    :param lM: Total Star Mass, :math:`\log(M_*)`
    :param mu: Mean Molecular Weight, :math:`\mu`
    :param m: Initial mass step from center of star, :math:`\log(dm)`
    :returns: Tuple of `(\log(r),\log(l),P,T)`, the four equation of state values at point m.
    """
    lrhoc = ldensity(logP=np.log10(Pc),logT=np.log10(Tc),mu=mu)
    rhoc = np.power(10,lrhoc)
    m = np.power(10,lm)
    epsilon = (epp(T=Tc,rho=rhoc,X=X,c=cfg) + eCNO(T=Tc,rho=rhoc,X=X,XCNO=XCNO,c=cfg))
    r = inner_radius(lrho=lrhoc,lm=lm)
    l = m * epsilon
    P = inner_pressure(Pc=Pc,rho=rhoc,m=m)
    T = inner_temperature(Tc=Tc,Pc=Pc,rho=rhoc,m=m,epsilon=epsilon,optable=optable,convective=convective)
    return (r,l,P,T)
    
    
def outer_boundary(R,L,M,mu,optable,Piter=False):
    r"""Load the outer boundary conditions for our star.
    
    :param R: Stellar Radius
    :param L: Total Luminosity
    :param M: Total Mass
    :param mu: Mean Molecular Weight :math:`\mu`
    :param function opfunc: A function which takes :math:`\log(T)` and :math:`\log(\rho)` and returns the rosseland mean opacity, :math:`\bar{\kappa}`
    :returns: Tuple of `(r,l,P,T)`, the four equation of state values at point m.
    
    """
    from scipy.optimize import fminbound
    
    r = R
    l = L
    T = outer_temperature(R=R,L=L)
    if Piter:
        if float(Piter) == 1.0:
            raise ValueError, "Can't iterate with multiplicative changes of =1.0"
        Pguess1 = find_pressure_guess(Pguess=1e1,T=T,mu=mu,optable=optable,incr=Piter)
        Pguess2 = find_pressure_guess(Pguess=1e10,T=T,mu=mu,optable=optable,incr=1/Piter)
    else:
        Pbound, Tbound = get_pressure_guess(T=T,mu=mu,optable=optable)
        Pguess2, Pguess1 = Pbound
    P, cost, ierr, numiter = fminbound(func=outer_pressure_cost,args=(R,T,M,mu,optable),x1=Pguess1,x2=Pguess2,full_output=True)
    if ierr:
        log.warning("No pressure convergance around P=[%g,%g]: Cost Remaining=%g" % (Pguess1,Pguess2,cost))
    return (r,l,P,T)
    
def log_outer_boundary(lR,lL,lM,mu,optable,Piter=False):
    r"""Load the outer boundary conditions for our star.
    
    :param lR: Stellar Radius, :math:`\log(R_*)`
    :param lL: Total Luminosity, :math:`\log(L_*)`
    :param lM: Total Mass, :math:`\log(M_*)`
    :param mu: Mean Molecular Weight :math:`\mu`
    :param function opfunc: A function which takes :math:`\log(T)` and :math:`\log(\rho)` and returns the rosseland mean opacity, :math:`\bar{\kappa}`
    :returns: Tuple of `(\log(r),\log(l),P,T)`, the four equation of state values at point m.
    
    """
    from scipy.optimize import fminbound
    
    lr = lR
    ll = L
    T = outer_temperature(R=np.power(10,lR),L=np.power(10,lL))
    if Piter:
        if float(Piter) == 1.0:
            raise ValueError, "Can't iterate with multiplicative changes of =1.0"
        Pguess1 = find_pressure_guess(Pguess=1e1,T=T,mu=mu,optable=optable,incr=Piter)
        Pguess2 = find_pressure_guess(Pguess=1e10,T=T,mu=mu,optable=optable,incr=1/Piter)
    else:
        Pbound, Tbound = get_pressure_guess(T=T,mu=mu,optable=optable)
        Pguess2, Pguess1 = Pbound
    P, cost, ierr, numiter = fminbound(func=outer_pressure_cost,args=(np.power(10,lR),T,np.power(10,lM),mu,optable),x1=Pguess1,x2=Pguess2,full_output=True)
    if ierr:
        log.warning("No pressure convergance around P=[%g,%g]: Cost Remaining=%g" % (Pguess1,Pguess2,cost))
    return (lr,ll,P,T)
    
    
    
def inner_radius(rho,m):
    r"""Inner radius to use, as a function of central density and initial step mass.
    
    .. math::
        
       r(m=m_0) = \frac{3}{4\pi \rho_c}^{1/3} m^{1/3}
    
    :param rho: Central density, :math:`\rho_c`
    :param m: Initial mass step from center of star.
    """
    return np.power(3/(4*np.pi*rho)*m,1/3)
    
def log_inner_radius(lrho,lm):
    r"""Inner radius to use, as a function of central density and initial step mass.
    
    .. math::
       
       \log(r) = \frac{1}{3}\log(\frac{3}{4\pi}) + \frac{1}{3}(m - \rho_c)
    
    :param rho: Central density, :math:`\rho_c`
    :param m: Initial mass step from center of star.
    """
    print lm - lrho
    return (1/3) * (np.log10(3/(4*np.pi)) + lm - lrho)
    
def inner_pressure(Pc,rho,m):
    r"""Inner pressure boundary of the system, given a specified density and mass step, and a central pressure.
    
    .. math::
        P(m=m_0) = P_c - \frac{3 G}{8 \pi} * \left(\frac{4 \pi}{3}\right)^{4/3} m^{2/3}
    
    :param Pc: Central Pressure, :math:`P_c`
    :param rho: Central density, :math:`\rho_c`
    :param m: Initial mass step from center of star.
    """
    from .constants import G
    Pdelt = ((3*G)/(8*np.pi) * np.power((4*np.pi)/3 * rho, 4/3) * np.power(m,2/3))
    return Pc - Pdelt
    
    
def inner_temperature(Tc,Pc,rho,m,epsilon,optable,convective=True):
    r"""Inner temperature for the core.
    
    .. math::
        
        \ln(T(m=m_0))-\ln(T_c) = - \left( \frac{\pi}{6} \right)^{1/3} \frac{\nabla_{AD} \rho_c^{4/3}}{P_c} m^{2/3}
        
    
    :param Pc: Central Pressure, :math:`P_c`
    :param Tc: Central Temperature, :math:`T_c`
    :param rho: Central density, :math:`\rho_c`
    :param m: Initial mass step from center of star.
    :param bool convective: Whether the star should be convective or not.
    
    """
    from .constants import gradT_ad, a, c
    if convective:
        lnT = np.log(Tc) - np.power(np.pi/6,1/3) * (gradT_ad * np.power(rho,4/3))/Pc * np.power(m,2/3)
        T = np.exp(lnT)
    else:
        optable.kappa(rho=rho,T=Tc)
        kappa = optable.retrieve()
        T4 = np.power(Tc,4) - 1/(2*a*c) * np.power((3/(4*np.pi)),2/3) * kappa * epsilon * np.power(rho,4/3) * np.power(m,2/3)
        T = np.power(T4,1/4)
    return T

def outer_pressure(R,M,kappa):
    r"""Outer pressure at the boundary condition.
    
    .. math::
    
        P(r=R) = \frac{GM}{R^2} \frac{2}{3} \frac{1}{\bar{\kappa}}
    
    :param R: Stellar Radius
    :param M: Total Mass
    :param kappa: Rosseland Mean Opacity
    
    """
    from .constants import G
    return (G * M)/np.power(R,2) * 2/3 * 1/kappa
    
def outer_pressure_cost(P,R,T,M,mu,optable):
    """Outer pressure cost function for the boundary condition.
    
    This function returns the difference between the pressure guessed and the pressure then calculated from the outer pressure equation implemented in :meth:`outer_pressure`
    
    """
    rho = density(P=P,T=T,mu=mu)
    optable.kappa(T=T,rho=rho)
    Pout = outer_pressure(R=R,M=M,kappa=optable.retrieve())
    return np.abs(Pout - P)
    
def outer_temperature(R,L):
    r"""Outer temperature at the boundary condition
    
    .. math::
        
        T = \left(\frac{L}{4\pi R^2 \sigma}\right)^{1/4}
        
    :param R: Stellar Radius
    :param L: Total luminosity
    """
    from .constants import sigmab
    return np.power(L/(4*np.pi*np.power(R,2)*sigmab), 1/4)
    
def get_pressure_guess(T,mu,optable):
    """Return an analytic estimate of the pressure boundary guesses, without iteration.
    
    This is imperfect, as it uses the corners of the table. The table isn't square, so it is still possible, for some choices of T and R that the requested opacity is off the charts.
    """
    from .constants import mh, kb, a
    optable.bounds()
    logRb,logTb = optable.retrieve()
    logTg = np.array([np.log10(T),np.log10(T)])
    pbounds = np.empty((2,2))
    optable.invert_points(logR=logRb,logT=logTg)
    logrho,logT = optable.retrieve().T
    Pbound = np.power(10,logrho) / (mu * mh) * kb * T + (a/3) * np.power(T,4)
    Tbound = np.power(10,logT)
    return np.vstack((Pbound,Tbound))
        
    
def find_pressure_guess(Pguess,T,mu,optable,incr=2):
    """Find a good initial pressure guess. Searches logarithmically, and checks to see if the guess is on the opacity table."""
    goodguess = False
    while not goodguess:
        rho = density(P=Pguess,T=T,mu=mu)
        points = np.array([np.log10(rho),np.log10(T)])
        optable.validate(points,Describe=True)
        goodguess,code,string = optable.retrieve()
        if not goodguess and code != 2**2:
            raise ValueError(string + " T=%g" % T)
        Pguess *= incr
    return Pguess
    
# For conformity with the "Numerical Recipies" names:

load1 = inner_boundary
load2 = outer_boundary