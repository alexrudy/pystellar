# -*- coding: utf-8 -*-
# 
#  newton.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-11-12.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 
from __future__ import division

import numpy as np
import scipy as sp
import scipy.linalg as la

class NRSolver(object):
    """Newton-Rhapson Solver"""
    def __init__(self, stars, config):
        super(NRSolver, self).__init__()
        self._stars = stars
        self._config = config
        
    
    @property
    def config(self):
        return self._config
        
    @property
    def stars(self):
        return self._stars
        
    def launch(self,nstar):
        """Launch integrations on nstar"""
        self.stars["Surface"][nstar].surface()
        self.stars["Center"][nstar].center()
        
    def fitgap(self,nstar):
        """Return the fit gap for nstar integrator threads"""
        fxs, fys, data = self.stars["Surface"][nstar].retrieve()
        fxs, fyc, data = self.stars["Center"][nstar].retrieve()
        return (fys[-1] - fyc[-1])
        
    def set_guesses(self,ys,nstar):
        """Set new guess values for nstar."""
        # Set new initial positions
        self.stars["Surface"][nstar].set_guesses(*ys)
        self.stars["Center"][nstar].set_guesses(*ys)
        self.stars["Surface"][nstar].release()
        self.stars["Center"][nstar].release()
        
    def fd_jacobian(self,y0):
        """Forward-difference jacobian from a y0 guess."""
        epsilon = float(self.config["System.NewtonRapson.dx"])
        dy = epsilon * y0
        nr = y0.size + 1
        y1 = np.empty((nr,y0.size))
        f = np.empty((nr,y0.size))
        jac = np.empty((y0.size,y0.size))
        y1[y0.size] = y0
        
        for i in range(y0.size):
            y1[i] = y0
            y1[i,i] = y0[i] + dy[i]
        
        for i in xrange(nr):
            self.set_guesses(y1[i],i)
            self.launch(i)
        
        for i in xrange(nr):
            f[i] = self.fitgap(i)
        print f
            
        for i in xrange(y0.size):
            jac[i] = (f[i] - f[-1]) #/dy
            
        return np.matrix(jac)
        
    def nrsolve(self):
        """Do the newton-rapson solution"""
        y0 = self._stars["master"][0].get_guesses()
        y0m = np.matrix(y0)
        print self.fd_jacobian(y0)
        # for n in xrange(self.config["System.NewtonRapson.niters"]):
            # jac = self.fd_jacobian(y0)
            # dy = - jac.I * y0m 