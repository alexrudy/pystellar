# -*- coding: utf-8 -*-
# 
#  constants.py
#  pystellar
#  Constants for a stellar model
#  
#  Created by Alexander Rudy on 2012-10-07.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 
"""
:mod:`constants` - Physical Constants
=====================================

All constants are in CGS units unless otherwise noted.

:mod:`astrophysics` constants
*****************************

.. data:: c
    
    The speed of light.
    
.. data:: mp
    
    The mass of a proton. It is assumed that :math:`m_p = m_h`.
    
.. data:: kb
    
    The Boltzmann Constant.
    
.. data:: G
    
    The gravitational constant.

Other constants
***************

"""

# cgs units!!
from astropysics.constants import c, mp, kb, G

sigmab = 5.6704e-5
"""Stephan-Boltzmann constant from http://en.wikipedia.org/wiki/Stefan–Boltzmann_constant"""

a = (4.0 * sigmab)/c 
"""radiation constant"""

mh = mp 
"""Mass of a hydrogen atom"""

gradT_ad = 0.4
"""Adiabatic Temperature Gradient"""


