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
        self.log.debug("Launching Star %d" % nstar)
        self.stars["Surface"][nstar].surface(plotting=plotting)
        self.stars["Center"][nstar].center(plotting=plotting)
        
    def fitgap(self,nstar):
        """Return the fit gap for nstar integrator threads"""
        fys, fxs, data = self.stars["Surface"][nstar].retrieve()
        fyc, fxs, data = self.stars["Center"][nstar].retrieve()
        fd = fys[-1] - fyc[-1]
        fn = np.min([fys[-1],fyc[-1]],axis=0)
        self.log.debug("Integration found R(fit)= %g - %g" % (fys[-1,0], fyc[-1,0]))
        if (np.abs(fd/fn) > float(self.config["System.NewtonRapson.maxfiterr"])).any():
            self.log.warning("Found fitting differences too large, declaring failure! %s / %s" % (fd,fn))
            fd = np.ones(fd.shape)*np.nan
        return fd
        
    def set_guesses(self,ys,nstar):
        """Set new guess values for nstar."""
        # Set new initial positions
        self.stars["Surface"][nstar].set_guesses(*ys)
        self.stars["Center"][nstar].set_guesses(*ys)
        self.stars["Surface"][nstar].release()
        self.stars["Center"][nstar].release()
        self.log.debug("Set Guesses for Star %d to %s" % (nstar,ys))
        
    def fd_jacobian(self,y0,f0,eps=None,depth=0):
        """Forward-difference jacobian from a y0 guess."""
        if eps is None:
            epsilon = float(self.config["System.NewtonRapson.Jac.dx"])
        else:
            epsilon = eps
        maxdepth = self.config["System.NewtonRapson.Jac.maxdepth"]
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
            self.log.debug("Jacobian Calcualtion Evaluated: %s" % ((f[i] - f0)/dy[i]))
            jac[i] = (f[i] - f0)/dy[i]
            
        if np.isfinite(jac).all():
            return np.matrix(jac.T),f
        elif depth < maxdepth:
            return self.fd_jacobian(y0,f0,eps=-0.1*epsilon,depth=depth+1)
        else:
            raise ValueError("Can't find convergent Jacobian")
        
    def nrsolve(self):
        """Do the newton-rapson solution"""
        # Reload previous initial guesses:
        if self.config.get("System.NewtonRapson.continue",False):
            y0 = np.load(self.config["System.NewtonRapson.file"])
            self.stars["master"][0].set_guesses(*y0)
            self.log.debug("Restarting Newton-Rapson Method from %s" % y0)
        else:
            y0 = self.stars["master"][0].get_guesses()
        
        
        tol = float(self.config["System.NewtonRapson.tol"])
        step_max = float(self.config["System.NewtonRapson.stepmax"])
        num_iter = self.config["System.NewtonRapson.niters"]
        converged = False
        
        # Initial values
        ny = y0.shape[0]
        self.set_guesses(y0,ny)
        self.launch(ny)
        f0 = self.fitgap(ny)
        
        # Tracking arrays
        fs = np.empty((num_iter,ny))
        ys = np.empty((num_iter,ny))
        cs = np.empty((num_iter,ny))
        ds = np.empty((num_iter,ny))
        
        if not np.isfinite(f0).all():
            self.log.critical("Can't converge on initial integration. Failing!")
            raise ValueError("Convergence failure in first integration")
        
        for n in xrange(num_iter):
            # Forward Difference Jacobian
            jac,fd1 = self.fd_jacobian(y0,f0)
            ff = 0.5 * np.dot(f0,f0)
            
            self.log.debug("Calculated Jacobian at R=%g, L=%g, Pc=%g, Tc=%g" % tuple(y0) )
            self.log.debug("Found Fitting Errors: %s" % f0)
            
            fs[n] = f0
            ys[n] = y0
            self.append_dashboard(np.array([n]),f0,line="fitting",marker='o',method='semilogy')
            self.append_dashboard(np.array([n]),y0,figure="guesses",line="guesses",marker='o')
            
            # Convergence Test
            c0 = np.abs(f0/y0)

            cs[n] = c0
            self.log.info("[%04d] Convergence Value: %.2g > %.2g" % (n,np.max(c0),tol))
            for i in range(ny):
                self.dashboard.append_data(np.array([n]),np.array([c0[i]]),figure='fitting',axes='ffpoint',
                    line='ffpoint-%d' % i,method='semilogy',marker='o')
            if np.max(c0) < tol:
                converged = True
                break
            
            # Matrix Math
            f0m = np.matrix(f0).T
            dy = -1*(jac.I * f0m).getA().T[0]
            dy *= float(self.config["System.NewtonRapson.stepeps"])
            
            ds[n] = dy
            self.append_dashboard(np.array([n]),dy,figure="adjustments",line="fitting-jac",marker="o")
            self.log.info("Full Step Size (dr,dl,dP,dT)=[%.4g, %.4g, %.4g, %.4g]" % tuple(dy))
            
            # Linear Search
            y0, f0 = self.linear_search(y0,f0,dy,ff)
            np.save(self.config["System.NewtonRapson.file"],y0)
            # Predictive update for next step
            self.append_dashboard(np.array([n]),f0,line="fitting-linear-search",marker='x',ms=6.0,method='semilogy')
            self.dashboard.update("fitting","adjustments","guesses")
            
        
        fs = fs[:n]
        ys = ys[:n]
        cs = cs[:n]
        ds = ds[:n]
        
        if converged:
            self.log.info("CONVERGENCE!! after %d iterations" % n)
        else:
            self.log.info("Reached Maximum Iterations: %d" % n)
        self.dashboard.update("fitting","adjustments","guesses")
        self.log.info("Guesses are: \nR=%g\nL=%g\nP=%g\nT=%g" % tuple(y0))
        # Save text files
        np.savetxt(self.config["System.Outputs.Data.FittingPoint"],fs)
        np.savetxt(self.config["System.Outputs.Data.Guesses"],ys)
        np.savetxt(self.config["System.Outputs.Data.Convergence"],cs)
        np.savetxt(self.config["System.Outputs.Data.Deltas"],ds)
        
        self.dashboard.save("fitting","fitting-points.pdf")
        self.dashboard.save("adjustments","guess-movement.pdf")
        self.dashboard.save("guesses","initial-conditions.pdf")
        self.dashboard.save("integration","final-integration.pdf")
        self.dashboard.save("integrationextras","final-integration-extra.pdf")
        
    
    def linear_search(self,x0,f0,dx,ff0):
        """Linear Search Method"""
        step_max = self.config["System.NewtonRapson.linearSearch.stepmax"] * max([np.sqrt(np.sum(np.power(x0,2))),x0.size])
        alpha = self.config["System.NewtonRapson.linearSearch.alpha"]
        step = np.sqrt(np.sum(dx**2))
        if step > step_max:
            self.log.debug("Reducing Step Size: %s" % dx)
            dx *= step_max/step
            self.log.warning("Reduced Step Size: %s" % dx)
        slope = np.sum(dx)
        if slope < 0:
            self.log.error("Roundoff Error in Search: dx=%s, slope=%.2g" % (dx,slope))
        test = np.max(np.abs(dx)/np.max(np.hstack((np.abs(x0),np.ones(x0.shape))),axis=0))
        alam_min = self.config["System.NewtonRapson.linearSearch.tolX"]/test
        alaminit = 1.0
        alam1 = alaminit
        converged = False
        convergence_steps = 0
        while not converged:
            convergence_steps += 1
            x1 = x0 + alam1 * dx
            self.log.debug("Moving guess by lambda %f" % alam1)
            self.set_guesses(x1,x1.size)
            self.launch(x1.size)
            if alam1 < alam_min:
                self.log.debug("Converged because lambda %f < lambda_min %f" % (alam1,alam_min))
                xr = x1
                converged = True
            f1 = self.fitgap(x1.size)
            ff1 = 0.5 * np.dot(f1,f1)
            if not np.isfinite(f1).all():
                alam0 = 0.9 * alam1
                self.log.warning("Non-finite integration failure in linear search")
                self.log.debug("Moving due to integration failure: %g -> %g" % (alam1,alam0))
            elif not converged and ff1 < ff0 + alpha * alam1 * slope:
                self.log.debug("Converged because ff %g < %g ff+alpha+lambda+slope (%g, %g, %g, %g)" % (ff1,ff0 + alpha * alam1 * slope,ff0,alpha,alam1,slope))
                xr = x1
                converged = True
            else:
                if alam1 == alaminit:
                    alam0 = - slope / (2.0 * (ff1 - ff0 * slope))
                    self.log.debug("Moving lambda on first step: %g -> %g" % (alam1,alam0))
                else:
                    rhs1 = ff1 - ff0 - alam1 * slope
                    rhs2 = ff2 - ff0 - alam2 * slope
                    a = (rhs1/(alam1**2) - rhs2/(alam2**2))/(alam1-alam2)
                    b = (-alam1 * rhs1/(alam1**2) + alam1 * rhs2/(alam2**2))/(alam1-alam2)
                    if a == 0.0:
                        alam0 = -slope/(2.0 * b)
                        self.log.debug("Moving only w/r/t b: %g -> %g" % (alam1,alam0))
                    else:
                        disc = b**2 - 3 * a * slope
                        if disc < 0.0:
                            alam0 = 0.5 * alam1
                            self.log.debug("Moving due to discriminant failure: %g -> %g" % (alam1,alam0))
                        elif b <= 0.0:
                            alam0 = (- b + np.sqrt(disc))/(3 * a)
                            self.log.debug("Moving due to positive root: %g -> %g" % (alam1,alam0))
                        else:
                            alam0 = -slope/(b + np.sqrt(disc))
                            self.log.debug("Moving due to negative root: %g -> %g" % (alam1,alam0))
                    if alam0 > 0.5 * alam1:
                        self.log.debug("Ensuring that we don't move by more than 0.5: %g -> %g" % (alam1,alam0))
                        alam0 = 0.5 * alam1
            if not converged:
                ff2 = ff1
                alam2 = alam1
                alam1 = 0.1 * alam1 if 0.1 * alam1 > alam0 else alam0
        self.log.info("Linear Search Converged after %d steps on %s" % (convergence_steps,xr))
        self.log.info("Final Step Size: %s" % (xr-x0))
        self.log.debug("Final Lambda Size: %g" % alam1)
        return xr,f1
    
    def append_dashboard(self,x,y,figure='fitting',line=None,**kwargs):
        """Append data to the dashboard"""
        if figure == 'fitting':
            y = np.abs(y)
        
        line = self.name if line is None else line
        for yi,name in zip(y,["radius","luminosity","pressure","temperature"]):
            yi = np.atleast_1d(yi)
            self.dashboard.append_data(x,yi,figure=figure,axes=name,line=line,**kwargs)
        