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
import logging

class NRSolver(object):
    """Newton-Rhapson Solver"""
    def __init__(self, stars, dashboard, config):
        super(NRSolver, self).__init__()
        self._stars = stars
        self._dashboard = dashboard
        self._config = config
        self.log = logging.getLogger(__name__)
        
    
    @property
    def dashboard(self):
        return self._dashboard
    
    @property
    def config(self):
        return self._config
        
    @property
    def stars(self):
        return self._stars
        
    def launch(self,nstar):
        """Launch integrations on nstar"""
        self.log.info("Launching Star %d" % nstar)
        self.stars["Surface"][nstar].surface()
        self.stars["Center"][nstar].center()
        
    def fitgap(self,nstar):
        """Return the fit gap for nstar integrator threads"""
        fys, fxs, data = self.stars["Surface"][nstar].retrieve()
        fyc, fxs, data = self.stars["Center"][nstar].retrieve()
        self.log.info("Integration found R(fit)= %g - %g" % (fys[-1,0], fyc[-1,0]))
        return (fys[-1] - fyc[-1])
        
    def set_guesses(self,ys,nstar):
        """Set new guess values for nstar."""
        # Set new initial positions
        self.stars["Surface"][nstar].set_guesses(*ys)
        self.stars["Center"][nstar].set_guesses(*ys)
        self.stars["Surface"][nstar].release()
        self.stars["Center"][nstar].release()
        self.log.info("Set Guesses for Star %d to %r" % (nstar,ys))
        
    def fd_jacobian(self,y0):
        """Forward-difference jacobian from a y0 guess."""
        
        epsilon = float(self.config["System.NewtonRapson.dx"])
        dy = epsilon * y0
        nr = y0.size + 1
        na = y0.size
        self.set_guesses(y0,na)
        self.launch(na)
        y1 = np.empty((na,na))
        f = np.empty((na,na))
        jac = np.empty((na,na))
        
        for i in range(na):
            y1[i] = y0
            y1[i,i] = y0[i] + dy[i]
        
        for i in xrange(na):
            self.set_guesses(y1[i],i)
            self.launch(i)
        
        for i in xrange(na):
            f[i] = self.fitgap(i)
        
        f0 = self.fitgap(na)
        
        for i in xrange(y0.size):
            jac[i] = (f[i] - f0)/dy
            
        return np.matrix(jac),f0
        
    def nrsolve(self):
        """Do the newton-rapson solution"""
        tol = self.config["System.NewtonRapson.tol"]
        y0 = self._stars["master"][0].get_guesses()
        
        for n in xrange(self.config["System.NewtonRapson.niters"]):
            y0m = np.matrix(y0).T
            jac,f0 = self.fd_jacobian(y0)
            ff = 0.5 * np.dot(f0,f0)
            self.log.info("Calculated Jacobian at R=%g, L=%g, Pc=%g, Tc=%g" % tuple(y0) )
            self.log.info("Found Fitting Errors: %r" % f0)
            self.append_dashboard(np.array([0]),f0,"fitting",marker='o')
            self.dashboard.update("fitting")
            if np.max(np.abs(f0)) < tol:
                break
            dy = (jac.I * y0m).getA().T[0]
            self.log.info("Moving fit by %r" % dy)
            y0 = self.linear_search(y0,f0,dy,ff)
            
        self.log.info("CONVERGENCE!! after %d iterations" % n)
        self.log.info("Guesses are: \nR=%g\nL=%g\nP=%g\nT=%g" % tuple(y0))
    
    def linear_search(self,x0,f0,dx,ff0):
        """docstring for linear_search"""
        step_max = self.config["System.NewtonRapson.linearSearch.stepmax"] * max([np.sqrt(np.sum(np.power(x0,2))),x0.size])
        alpha = self.config["System.NewtonRapson.linearSearch.alpha"]
        step = np.sqrt(np.sum(np.power(-dx,2)))
        if step > step_max:
            dx *= step_max/step
        slope = np.sum(-dx * f0)
        if slope > 0:
            self.log.error("Roundoff Error in Search: dx=%r f0=%r" % (dx,f0))
        test = np.max(np.abs(dx)/np.max(np.hstack((np.abs(x0),np.ones(x0.shape))),axis=0))
        alam_min = self.config["System.NewtonRapson.linearSearch.tolX"]/test
        alam1 = 1.0
        converged = False
        convergence_steps = 0
        while not converged:
            convergence_steps += 1
            x1 = x0 + alam1 * dx
            self.set_guesses(x1,x1.size)
            assert False
            self.launch(x1.size)
            if alam1 < alam_min:
                xr = x0
                converged = True
            f1 = self.fitgap(x1.size)
            ff1 = 0.5 * np.dot(f1,f1)
            if not converged and ff1 < ff0 + alpha * alam1 * slope:
                xr = x0
                converged = True
            else:
                if alam1 == 1.0:
                    alam0 = - slope / (2.0 * (ff1 - ff0 * slope))
                else:
                    rhs1 = ff1 - ff0 - alam1 * slope
                    rhs2 = ff2 - ff0 - alam2 * slope
                    a = (rhs1/(alam1**2) - rhs2/(alam2**2))/(alam1-alam2)
                    b = (-alam1 * rhs1/(alam1**2) + alam1 * rhs2/(alam2**2))/(alam1-alam2)
                    if a == 0.0:
                        alam0 = -slope/(2.0 * b)
                    else:
                        disc = b**2 - 3 * a * slope
                        if disc < 0.0:
                            alam0 = 0.5 * alam
                        elif b <= 0.0:
                            alam0 = (- b + np.sqrt(disc))/(3 * a)
                        else:
                            alam0 = -slope/(b + np.sqrt(disc))
                    if alam0 > 0.5 * alam1:
                        alam0 = 0.5 * alam1
            ff2 = ff1
            alam2 = alam1
            alam1 = 0.1 * alam1 if 0.1 * alam1 > alam0 else alam0
        self.log.debug("Linear Search Converged after %d steps on %g" % (convergence_steps,xr))
        return xr
    
    def append_dashboard(self,x,y,line=None,**kwargs):
        """Append data to the dashboard"""
        line = self.name if line is None else line
        for yi,name in zip(y,["radius","luminosity","pressure","temperature"]):
            self.dashboard.append_data(x,yi,figure="fitting",axes=name,line=line,**kwargs)
        