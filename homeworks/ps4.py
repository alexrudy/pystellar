#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  energy_plot.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-11-12.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

import numpy as np
import matplotlib.pyplot as plt

from pystellar.energy import epp, eCNO

from pystellar.opacity import OpacityTable
from pystellar.threading import ObjectThread

from pyshell.config import DottedConfiguration

Config = DottedConfiguration()
Config.dn = DottedConfiguration
Config.load("pystellar/Star.yml")
Config.load("Star.yml")
X = 0.700
Y = 0.280
XCNO = 0.7 * (1 - X - Y)
rho = 80

T = np.logspace(6.8,7.8,100)

epsilon_pp = epp(T=T,rho=rho,X=X,c=Config["Data.Energy"])
epsilon_CNO = eCNO(T=T,rho=rho,X=X,XCNO=XCNO,c=Config["Data.Energy"])

T1 = 18 * np.power(10,6)
epsilon_pp1 = epp(T=T1,rho=rho,X=X,c=Config["Data.Energy"])
epsilon_CNO1 = eCNO(T=T1,rho=rho,X=X,XCNO=XCNO,c=Config["Data.Energy"])

print u"X        = %g" % X
print u"XCNO     = %g" % XCNO
print u"ρ        = %g" % rho
print u"T        = %g" % T1
print u"εCNO     = %g" % epsilon_CNO1
print u"εpp      = %g" % epsilon_pp1
print u"ε        = %g" % (epsilon_pp1 + epsilon_CNO1)
print u"εpp/εCNO = %g" % (epsilon_pp1/epsilon_CNO1)


plt.loglog(T,epsilon_CNO+epsilon_pp,'-',label=r"$\epsilon$")
plt.loglog(T,epsilon_pp,'-.',lw=4.0,label=r"$\epsilon_\textrm{pp}$")
plt.loglog(T,epsilon_CNO,'--',lw=4.0,label=r"$\epsilon_\textrm{CNO}$")
plt.xlim((10**6.8,10**7.8))
plt.ylim((10**-3,10**6))
plt.xlabel(r"$T$")
plt.ylabel(r"$\epsilon$")
plt.legend()
plt.savefig("epsilon.pdf")

