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

from .density import density

def load_inner(Pc,Tc,M,mu,m):
    r"""Load our inner boundary conditions for our star.
    
    :param Pc: Central Pressure, :math:`P_c`
    :param Tc: Central Temperature, :math:`T_c`
    :param M: Total Star Mass
    :param mu: Mean Molecular Weight, :math:`\mu`
    :param m: Initial mass step from center of star.
    :returns: Tuple of `(r,l,P,T)`, the four equation of state values at point m.
    
    .. todo::
        I'm not sure if this is right... the local luminosity should probably depend on nuclear energy generation rates.
    
    """
    rhoc = density(P=Pc,T=Tc,mu=mu)
    r = inner_radius(rho=rhoc,m=m)
    l = 0 #TODO: I'm not sure about this!
    P = inner_pressure(Pc=Pc,rho=rhoc,m=m)
    T = inner_temperature(Tc=Tc,Pc=Pc,rho=rho,m=m)
    return (r,l,P,T)
    
def load_outer(R,L,M,mu,opfunc):
    r"""Load the outer boundary conditions for our star.
    
    :param R: Stellar Radius
    :param L: Total Luminosity
    :param M: Total Mass
    :param mu: Mean Molecular Weight :math:`\mu`
    :param function opfunc: A function which takes :math:`\log(T)` and :math:`\log(\rho)` and returns the rosseland mean opacity, :math:`\bar{\kappa}`
    :returns: Tuple of `(r,l,P,T)`, the four equation of state values at point m.
    
    .. todo::
        Fix recursive dependency between opacity and pressure for initial conditions. Maybe provide an initial rosseland mean opacity?
    
    """
    r = R
    l = L
    Phuh = 0 #TODO: We need some pressure here as our initial guess, in order to get density, in order to get pressure.
    T = outer_temperature(R=R,L=L)
    rho = density(P=Phuh,T=T,mu=mu) #TODO: Finish this function...
    kappa = opfunc(logT=np.log(T),logrho=np.log(rho))
    P = outer_pressure(R=R,M=M,kappa=kappa)
    return (r,l,P,T)
    
def inner_radius(rho,m):
    r"""Inner radius to use, as a function of central density and initial step mass.
    
    .. math::
        
       r(m=m_0) = \frac{3}{4\pi \rho_c}^{1/3} m^{1/3}
    
    :param rho: Central density, :math:`\rho_c`
    :param m: Initial mass step from center of star.
    """
    return np.power(3/(4*np.pi*rho)*m,1/3)
    
def inner_pressure(Pc,rho,m):
    r"""Inner pressure boundary of the system, given a specified density and mass step, and a central pressure.
    
    .. math::
        P(m=m_0) = P_c - \frac{3 G}{8 \pi} * \left(\frac{4 \pi}{3})\right)^{4/3} m^{2/3}
    
    :param Pc: Central Pressure, :math:`P_c`
    :param rho: Central density, :math:`\rho_c`
    :param m: Initial mass step from center of star.
    """
    from .constants import G
    return Pc-(3*G)/(8*np.pi) * np.power((4*np.pi)/3 * rho, 4/3) * np.power(m,2/3)
    
def inner_temperature(Tc,Pc,rho,m,convective=True):
    r"""Inner temperature for the core.
    
    .. math::
        
        \ln(T(m=m_0))-\ln(T_c) = - \left( \frac{\pi}{6} \right)^{1/3} \frac{\nabla_{AD} \rho_c^{4/3}}{P_c} m^{2/3}
        
    
    :param Pc: Central Pressure, :math:`P_c`
    :param Tc: Central Temperature, :math:`T_c`
    :param rho: Central density, :math:`\rho_c`
    :param m: Initial mass step from center of star.
    :param bool convective: Whether the star should be convective or not.
    
    """
    if not convective:
        raise NotImplementedError("Cannot calculate radiative core stars right now.")
    D_AD = 0.2 #Adiabatic Temperature Gradient
    lnT = np.log(Tc) - np.power(np.pi/6,1/3) * (D_AD* np.power(rho,4/3))/Pc * np.power(m,2/3)
    return np.exp(lnT)

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
    
def outer_temperature(R,L):
    r"""Outer temperature at the boundary condition
    
    .. math::
        
        T = \left(\frac{L}{4\pi R^2 \sigma}\right)^{1/4}
        
    :param R: Stellar Radius
    :param L: Total luminosity
    """
    return np.power(L/(4*np.pi*np.power(R,2)), 1/4)
    
# For conformity with the "Numerical Recipies" names:

load1 = load_inner
load2 = load_outer