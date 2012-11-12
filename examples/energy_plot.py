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
Config.load("Star.yml")
X = 1.000
Y = 0.000
rho = 1
Opacity = ObjectThread(OpacityTable,ikwargs=dict(fkey='GN93hz',X=X,Y=Y,snap=True),locking=True,timeout=None)
Opacity.start()

T = np.logspace(6.5,7.5,100)

epsilon_pp = epp(T=T,rho=rho,X=X,c=Config["Data.Energy"])
epsilon_CNO = eCNO(T=T,rho=rho,X=X,XCNO=0.01,c=Config["Data.Energy"])

plt.loglog(T,epsilon_CNO+epsilon_pp,'-',label='Total')
plt.loglog(T,epsilon_pp,':',label="PP")
plt.loglog(T,epsilon_CNO,'--',label="CNO")
plt.legend()
plt.show()

Opacity.stop()
