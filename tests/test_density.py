# -*- coding: utf-8 -*-
# 
#  test_density.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-07.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 
from __future__ import division

from pystellar.density import ldensity, mmw, lbeta

Xs    = [0,     0.70]
Ys    = [0.98,  0.28]
logTs = [7.55,  6.91]
logPs = [16.85,16.87]


class test_density(object):
    """Desntiy Calculation"""
    
    def test_solar_mu(self):
        """Solar values work for mu(X,Y): X=0.70, Y=0.28, log(T)=6.91, log(P)=16.87"""
        X = 0.70
        Y = 0.28
        logT = 6.91
        logP = 16.87
        mu = mmw(X=X,Y=Y)
        mu_sol = 0.617283950617284
        assert mu - mu_sol < 1e-10, u"μ mismatch: %g, %g" % (mu,mu_sol)
        
    def test_solar_beta(self):
        """Solar values work for beta(log(T),log(P),mu) : X=0.70, Y=0.28, log(T)=6.91, log(P)=16.87"""
        X = 0.70
        Y = 0.28
        logT = 6.91
        logP = 16.87
        beta = lbeta(logT=logT,logP=logP,X=X,Y=Y)
        beta_sol = 0.99985149819598074
        assert beta - beta_sol < 1e-10, u"β mimatch: %g, %g" % (beta, beta_sol)
        
        
    def test_solar_rho(self):
        """Solar values work for rho(log(T),log(P),mu) : X=0.70, Y=0.28, log(T)=6.91, log(P)=16.87"""
        X = 0.70
        Y = 0.28
        logT = 6.91
        logP = 16.87
        dens = ldensity(logT=logT,logP=logP,X=X,Y=Y)
        dens_sol = 68.189717584403212
        assert dens - dens_sol < 1e-10, u"ρ mismatch: %g, %g" % (dens, dens_sol)