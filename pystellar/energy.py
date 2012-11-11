# -*- coding: utf-8 -*-
# 
#  energy.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-11-10.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

import numpy as np

def epp(T,rho,X,c):
    """docstring for epp"""
    T6 = T / 1e6
    psi = 1.0
    g11 = 1 + float(c["g11.A"]) * np.power(T6,1/3) + float(c["g11.B"]) * np.power(T6,2/3) + float(c["g11.C"])*T6
    epsilon_pp = float(c["epp.A"]) * psi * float(c["epp.f11"]) * g11 * rho * np.power(X,2) * np.power(T6,-2/3) * np.exp(-float(c["epp.B"])/np.power(T6,1/3))
    return epsilon_pp
    
def eCNO(T,rho,X,XCNO,c):
    """docstring for fname"""
    T6 = T / 1e6
    g141 = 1 + float(c["g141.A"]) * np.power(T6,1/3) - float(c["g141.B"]) * np.power(T6,2/3) - float(c["g141.C"]) * T6
    epsilon_CNO = float(c["eCNO.A"])  * g141 * rho * X * XCNO * np.power(T6,-2/3) * np.exp(-float(c["eCNO.B"])/np.power(T6,1/3))
    return epsilon_CNO
    