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
        
    def launch(self,nstar,plotting=True):
        """Launch integrations on nstar"""
        self.log.info("Launching Star %d" % nstar)
        self.stars["Surface"][nstar].surface(plotting=plotting)
        self.stars["Center"][nstar].center(plotting=plotting)
        
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
        self.log.debug("Set Guesses for Star %d to %r" % (nstar,ys))
        
    def fd_jacobian(self,y0,f0):
        """Forward-difference jacobian from a y0 guess."""
        
        epsilon = float(self.config["System.NewtonRapson.dx"])
        dy = epsilon * y0
        nr = y0.size + 1
        na = y0.size
        y1 = np.empty((na,na))
        f = np.empty((na,na))
        jac = np.empty((na,na))
        
        for i in range(na):
            y1[i] = y0
            y1[i,i] = y0[i] + dy[i]
        
        for i in xrange(na):
            self.set_guesses(y1[i],i)
            self.launch(i,plotting=False)
        
        for i in xrange(na):
            f[i] = self.fitgap(i)
        
        for i in xrange(y0.size):
            self.log.info("Jacobian Calcualtion Evaluated: %s" % ((f[i] - f0)/dy[i]))
            jac[i] = (f[i] - f0)/dy[i]
            
        return np.matrix(jac.T),f
        
    def nrsolve(self):
        """Do the newton-rapson solution"""
        tol = float(self.config["System.NewtonRapson.tol"])
        step_max = float(self.config["System.NewtonRapson.stepmax"])
        y0 = self._stars["master"][0].get_guesses()
        converged = False
        self.set_guesses(y0,y0.size)
        self.launch(y0.size)
        f0 = self.fitgap(y0.size)
        for n in xrange(self.config["System.NewtonRapson.niters"]):
            y0m = np.matrix(y0).T
            jac,fd1 = self.fd_jacobian(y0,f0)
            ff = 0.5 * np.dot(f0,f0)
            self.log.info("Calculated Jacobian at R=%g, L=%g, Pc=%g, Tc=%g" % tuple(y0) )
            self.log.info("Found Fitting Errors: %s" % f0)
            self.append_dashboard(np.array([n]),f0,line="fitting",marker='o')
            self.append_dashboard(np.array([n]),y0,figure="guesses",line="guesses",marker='o')
            if np.max(np.abs(f0)) < tol:
                converged = True
                break
            dy = -1*(jac.I * y0m).getA().T[0]
            dy *= float(self.config["System.NewtonRapson.stepeps"])
            self.log.info("Step Size Required: %s" % dy)
            while (np.abs(dy) > y0 * step_max).any():
                self.log.warning("Step size is too large! %s" % dy)
                dy *= step_max/np.sqrt(np.sum(dy**2))
            self.append_dashboard(np.array([n]),dy,figure="adjustments",line="fitting-jac",marker="o")
            self.dashboard.update("fitting","adjustments","guesses")
            self.log.info("Moving fit by total dy= %g, %g, %g, %g" % tuple(dy))
            y1, f1 = self.linear_search(y0,f0,dy,ff)
            self.append_dashboard(np.array([n]),f1,line="fitting-linear-search",marker='x',ms=6.0)
            f0 = f1
            y0 = y1
            self.dashboard.update("fitting")
            
        if converged:
            self.log.info("CONVERGENCE!! after %d iterations" % n)
        else:
            self.log.info("Reached Maximum Iterations: %d" % n)
        self.log.info("Guesses are: \nR=%g\nL=%g\nP=%g\nT=%g" % tuple(y0))
    
    def linear_search(self,x0,f0,dx,ff0):
        """Linear Search Method"""
        step_max = self.config["System.NewtonRapson.linearSearch.stepmax"] * max([np.sqrt(np.sum(np.power(x0,2))),x0.size])
        alpha = self.config["System.NewtonRapson.linearSearch.alpha"]
        step = np.sqrt(np.sum(dx**2))
        if step > step_max:
            self.log.warning("Linear Search Step size is too large! %r" % dy)
            dx *= step_max/step
        slope = np.sum(-dx * f0)
        if slope > 0:
            self.log.error("Roundoff Error in Search: dx=%r f0=%r" % (dx,f0))
        test = np.max(np.abs(dx)/np.max(np.hstack((np.abs(x0),np.ones(x0.shape))),axis=0))
        alam_min = self.config["System.NewtonRapson.linearSearch.tolX"]/test
        alaminit = 0.6
        alam1 = alaminit
        converged = False
        convergence_steps = 0
        while not converged:
            convergence_steps += 1
            x1 = x0 + alam1 * dx
            self.set_guesses(x1,x1.size)
            self.launch(x1.size)
            if alam1 < alam_min:
                xr = x1
                converged = True
            f1 = self.fitgap(x1.size)
            ff1 = 0.5 * np.dot(f1,f1)
            if not converged and ff1 < ff0 + alpha * alam1 * slope:
                xr = x1
                converged = True
            else:
                if alam1 == alaminit:
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
                            alam0 = 0.5 * alam1
                        elif b <= 0.0:
                            alam0 = (- b + np.sqrt(disc))/(3 * a)
                        else:
                            alam0 = -slope/(b + np.sqrt(disc))
                    if alam0 > 0.5 * alam1:
                        alam0 = 0.5 * alam1
            if not converged:
                ff2 = ff1
                alam2 = alam1
                alam1 = 0.1 * alam1 if 0.1 * alam1 > alam0 else alam0
        self.log.info("Linear Search Converged after %d steps on %s" % (convergence_steps,xr))
        self.log.info("Final Step Size: %s" % (xr-x0))
        return xr,f1
    
    def append_dashboard(self,x,y,figure='fitting',line=None,**kwargs):
        """Append data to the dashboard"""
        line = self.name if line is None else line
        for yi,name in zip(y,["radius","luminosity","pressure","temperature"]):
            self.dashboard.append_data(x,yi,figure=figure,axes=name,line=line,**kwargs)
        