# -*- coding: utf-8 -*-
# 
#  star.py
#  pystellar
#  
#  Created by Jaberwocky on 2012-11-09.
#  Copyright 2012 Jaberwocky. All rights reserved.
# 
from __future__ import division

import numpy as np
import scipy as sp

from AstroObject.config import DottedConfiguration

import logging
from multiprocessing import Process
from StringIO import StringIO

from .opacity import OpacityTable
from .threading import ObjectThread, ObjectManager
from .initial import inner_boundary, log_inner_boundary, log_outer_boundary, outer_boundary
from .stellar import derivatives, log_derivatives, radiative_gradient, grad, dldm
from .density import mmw, density
from .integrator import integrate

class Star(object):
    """A simple configured star object."""
    def __init__(self, config=None, optable_args=None, dashboard_args=None):
        super(Star, self).__init__()
        self.log = logging.getLogger(__name__)
        self.data_log = logging.getLogger("data")
        self.telem_log = logging.getLogger("telemetry")
        self.telem_log.propagate = False
        self._call_count = 0
                
        self._config = DottedConfiguration(config)
        self._config.dn = DottedConfiguration
        
        self._integrator_name = self.config.get("System.Integrator.Method","scipy")
        self._integrator = getattr(self,self._integrator_name)
        self._logmode = self.config.get("System.Integrator.LogScale",True)
        np.set_printoptions(**self.config["System.Numpy.PrintOptions"])
        
        if optable_args is None:
            self.log.debug("Starting Opacity Table from Scratch")
            self._opacity = ObjectThread(OpacityTable,
                ikwargs=dict(fkey=self._config["Data.Opacity.Config"],X=self.X,Y=self.Y,
                    snap=self._config["System.Opacity.Snap"],error=self._config["System.Opacity.Error"]),
                locking=True,timeout=self._config["System.Opacity.Timeout"])
            self.opacity.start()
            self.log.debug("Started Opacity Table From Scratch")
        else:
            self.log.debug("Starting Opacity Table from Arguments")
            self._opacity = ObjectManager(**optable_args)
            self.log.debug("Started Opacity Table from Arguments")
        
        if dashboard_args is not None:
            self._dashboard = ObjectManager(**dashboard_args)
        
        from astropysics.constants import Rs, Lsun, Ms
        self.Pc_Guess = float(self._config["Star.Initial.Pc"])
        self.Tc_Guess = float(self._config["Star.Initial.Tc"])
        self.R_Guess = float(self._config["Star.Initial.Rs"]) * Rs
        self.L_Guess = float(self._config["Star.Initial.Ls"]) * Lsun
        self.dm_Guess = float(self._config["Star.Initial.dm"]) * Ms
        self.fp = float(self._config["Star.Integration.fp"]) * self.M
        
        self._update_frequency = self.config["System.Dashboard.Update"]
        self._save_frequency = self.config["System.Outputs.Update"]
        
        
    def set_fitting_point(self,point):
        """Set a new fitting point for this routine"""
        self.fp = point
        
    def set_guesses(self,R,L,Pc,Tc):
        """docstring for set_guesses"""
        self.Pc_Guess = Pc
        self.Tc_Guess = Tc
        self.R_Guess = R
        self.L_Guess = L
        
    def integral(self,y,x,i):
        """docstring for integral"""
        x = np.asarray(x)
        y = np.asarray(y)
        self._call_count += 1            
        
        dy = self.fprime(xs=x,ys=y,mu=self.mu,optable=self.opacity,X=self.X,XCNO=self.Z,cfg=self.config["Data.Energy"])[:,0]
        
        if self._logmode:
            x = np.power(10,x)
        update = self._update_frequency != 0 and self._call_count % self._update_frequency == 0
        telem  = self._save_frequency != 0 and self._call_count % self._save_frequency == 0
        if (y < 0).any() or update or telem:
            rho = density(P=y[2],T=y[3],mu=self.mu)
            eps = dy[1]
            self.opacity.kappa(T=y[3],rho=rho)
            rgrad = radiative_gradient(T=y[3],P=y[2],l=y[1],m=x,rho=rho,optable=self.opacity)
            agrad = grad(rgrad)
            self.opacity.kappa(T=y[3],rho=rho)
            kappa = self.opacity.retrieve()
        if (y < 0).any():
            self.log.warning(u"%s y<0 at: \nx=%r, \ny=%r, \ndy=%r, \nρ=%g \nε=%g \n∇=%g \nκ=%g" % (i,x,y,dy,rho,eps,agrad,kappa[0]))
        if update:
            self.dashboard.insert_data(x,y,i)
            self.dashboard.add_density(x,rho,i)
            self.dashboard.add_epsilon(x,dy[1],i)
            self.log.info(u"%s %d calls at: \nx=%r, \ny=%r, \ndy=%r, \nρ=%g \nε=%g \n∇=%g \nκ=%g" % (i,self._call_count,x,y,dy,rho,eps,agrad,kappa[0]))
            self.dashboard.update()
        if telem:
            self.telem_log.info(u"%r %r %r %g %g %g %g" % (x,y,dy,rho,eps,agrad,kappa[0]))
        return dy
        
    @property
    def integrator(self):
        """The actual integrator function!"""
        return self._integrator
        
    @property
    def fprime(self):
        """The correct fprime function for integration"""
        if self._logmode:
            return log_derivatives
        else:
            return derivatives
        
        
    def show_surface_start(self):
        """Show a detailed view of the start of integration at the center"""
        self.fp = (1-1e-3) * self.M
        self._config["System.Outputs.Size"] = 10
        return self.surface()
        
    def show_center_start(self):
        """Show a detailed view of the start of integration at the center"""
        self.fp = 1e30
        self._config["System.Outputs.Size"] = 10
        return self.center()
        
    def center(self):
        """Run the center integration."""
        self.log.debug("Getting Inner Boundaries")
        center_ic = inner_boundary(
            Pc=self.Pc_Guess,Tc=self.Tc_Guess,M=self.M,mu=self.mu,m=self.dm_Guess,
            optable=self.opacity,X=self.X,XCNO=self.XCNO,cfg=self.config["Data.Energy"],convective=self.config["Star.Initial.Convective"])
        if self._logmode:
            integrator = "LogCenter"
            ms = np.linspace(np.log10(self.dm_Guess),np.log10(self.fp),self._config["System.Outputs.Size"])
        else:
            integrator = "Center"
            ms = np.logspace(np.log10(self.dm_Guess),np.log10(self.fp),self._config["System.Outputs.Size"])
            
        self.log.debug("Inner Conditions (%s): x=%g, y=%r" % (integrator,self.dm_Guess,np.array(center_ic)))
        self.log.debug("Starting %s Integration" % integrator)
        
        return self.integrator(ms,center_ic,integrator)
        
        
    def surface(self):
        """Run an integration from the surface to the inner edge"""
        self.log.debug("Getting Inner Boundaries")
        outer_ic = outer_boundary(R=self.R_Guess,L=self.L_Guess,M=self.M,mu=self.mu,optable=self.opacity)
        (r,l,P,T) = outer_ic
        
        if self._logmode:
            integrator = "LogSurface"
            ms = np.linspace(np.log10(self.fp),np.log10(self.M),self._config["System.Outputs.Size"])[::-1]
        else:
            integrator = "Surface"
            ms = np.logspace(np.log10(self.fp),np.log10(self.M),self._config["System.Outputs.Size"])[::-1]
        
        self.log.debug("Surface Conditions (%s): x=%g, y=%r" % (integrator,self.M,np.array(outer_ic)))
        self.log.debug("Starting %s Integration" % integrator)
        return self.integrator(ms,outer_ic,integrator)
        
    
    def pystellar(self,xs,ics,integrator):
        """Run an integration from the central point to the outer edge."""
        xs,ys,xc,yc = integrate(self.integral,xs,ics,args=(integrator,),**self.config["System.Integrator.PyStellar"][integrator]["Arguments"])
        self.log.debug("Finished %s Integration" % integrator)
        return ys, xs, None
        
    def scipy(self,xs,ics,integrator):
        """Run an integration from the central point to the outer edge."""
        self.log.debug("Calling %s Scipy Integration" % integrator)
        import scipy.integrate
        ys, data = scipy.integrate.odeint(self.integral,ics,xs,args=(integrator,),
            full_output=True,**self.config["System.Integrator.Scipy"][integrator]["Arguments"])
        self.log.debug("Finished %s Integration: %d derivative calls." % (integrator,self._call_count))
        if self._logmode:
            xs = np.power(10,xs)
        rho = density(P=ys[:,2],T=ys[:,3],mu=self.mu)
        eps = dldm(T=ys[:,3],rho=rho,X=self.X,XCNO=self.XCNO,cfg=self.config["Data.Energy"])
        self.opacity.kappa(T=ys[:,3],rho=rho)
        rgrad = radiative_gradient(T=ys[:,3],P=ys[:,2],l=ys[:,1],m=xs,rho=rho,optable=self.opacity)
        agrad = grad(rgrad)
        self.opacity.kappa(T=ys[:,3],rho=rho)
        kappa = self.opacity.retrieve()
        all_data = np.vstack(map(np.atleast_2d,(xs,ys.T,rho,eps,rgrad,agrad,kappa))).T
        stream = StringIO()
        np.savetxt(stream,all_data)
        self.data_log.info(stream.getvalue())
        stream.close()
        self.dashboard.replace_data(xs,ys,integrator)

        # self.opacity.kappa(T=ys[:,3],rho=rho)
        # eps = radiative_gradient(T=ys[:,3],P=ys[:,2],l=ys[:,1],m=xs,rho=rho,optable=self.opacity)
        self.dashboard.add_density(xs,rho,integrator,append=False)
        self.dashboard.add_epsilon(xs,eps,integrator,append=False)
        self.dashboard.update()
        self.log.debug("Plotted %s Integration" % integrator)
        return ys, xs, data
    
    
    @property
    def dashboard(self):
        """Dashboard Thread"""
        return self._dashboard
    
    
    @property
    def opacity(self):
        """Opacity table"""
        return self._opacity
        
    # Compositon Accessor Functions
    @property
    def X(self):
        """Fraction of Hydrogen"""
        return self._config["Star.Composition.X"]
        
    @property
    def Y(self):
        """Fraction of Helium"""
        return self._config["Star.Composition.Y"]
        
    @property
    def dXc(self):
        """Excess fraction of Carbon"""
        return self._config["Star.Composition.dXc"]
        
    @property
    def dXo(self):
        """Excess fraction of Oxygen"""
        return self._config["Star.Composition.dXo"]
        
    @property
    def Z(self):
        """Fraction of Heavy Elements"""
        return 1.0 - (self.X + self.Y + self.dXc + self.dXo)
    
    @property
    def XCNO(self):
        """Fraction of Carbon, Oxygen and Nitrogen"""
        return self.config["Star.Composition.ZdXCNO"] * self.Z + self.dXc + self.dXo
    
    @property
    def mu(self):
        """Mean molecular weight of this star."""
        return mmw(X=self.X,Y=self.Y)
    
    @property
    def M(self):
        """Total Mass of the Star"""
        from astropysics.constants import Ms
        return self._config["Star.Properties.M"] * Ms
    
    # Configuration Functions
    @property
    def config(self):
        """The configuration for this star."""
        return self._config
        
    def stop(self):
        """End this star!"""
        pass
        