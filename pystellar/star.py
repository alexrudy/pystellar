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

from pystellar.opacity import OpacityTable
from pystellar.threading import ObjectThread, ObjectManager

from pystellar.initial import inner_boundary, outer_boundary
from pystellar.stellar import derivatives
from pystellar.density import mmw

class Star(object):
    """A simple configured star object."""
    def __init__(self, filename="Star.yml", optable_args=None):
        super(Star, self).__init__()
        self.log = logging.getLogger(__name__)
        self._filename = filename
        self._config = DottedConfiguration()
        self._config.load(self._filename)
        self._config.dn = DottedConfiguration
        self.log.debug("Optable Args: %r" % optable_args)
        if optable_args is None:
            self.log.debug("Starting Opacity Table from Scratch")
            self._opacity = ObjectThread(OpacityTable,ikwargs=dict(fkey=self._config["Data.Opacity.Config"],X=self.X,Y=self.Y,snap=True),locking=True,timeout=None)
            self.opacity.start()
            self.log.debug("Started Opacity Table From Scratch")
        else:
            self.log.debug("Starting Opacity Table from Arguments")
            self._opacity = ObjectManager(**optable_args)
            self.log.debug("Started Opacity Table from Arguments")
        
        from astropysics.constants import Rs, Lsun, Ms
        self.Pc_Guess = float(self._config["Star.Initial.Pc"])
        self.Tc_Guess = float(self._config["Star.Initial.Tc"])
        self.R_Guess = float(self._config["Star.Initial.Rs"]) * Rs
        self.L_Guess = float(self._config["Star.Initial.Ls"]) * Lsun
        self.dm_Guess = float(self._config["Star.Initial.dm"]) * Ms
        self.dummy_epsilon = 100
        self.fp = 1e18
        
    def set_fitting_point(self,point):
        """Set a new fitting point for this routine"""
        self.fp = point
        
    def set_guesses(self,R,L,Pc,Tc):
        """docstring for set_guesses"""
        self.Pc_Guess = Pc
        self.Tc_Guess = Tc
        self.R_Guess = R
        self.L_Guess = L
        
    def kill(self):
        """Kill any processes we own!"""
        if isinstance(self._opacity,Process):
            self._opacity.kill()
            
    def stop(self):
        """Stop any processes we own!"""
        if isinstance(self._opacity,Process):
            self._opacity.stop()
        
    def integral(self,y,x):
        """docstring for integral"""
        return derivatives(xs=x,ys=y,mu=self.mu,epsilon=self.dummy_epsilon,optable=self.opacity)[0,:]
        
    def center(self):
        """Run an integration from the central point to the outer edge."""
        import scipy.integrate
        self.log.debug("Getting Inner Boundaries")
        center_ic = inner_boundary(Pc=self.Pc_Guess,Tc=self.Tc_Guess,M=self.M,mu=self.mu,m=self.dm_Guess,optable=self.opacity,epsilon=self.dummy_epsilon)
        ms = np.logspace(np.log10(self.dm_Guess),np.log10(self.fp),500)
        self.log.debug("Starting Integration")
        ys, data = scipy.integrate.odeint(self.integral,center_ic,ms,full_output=True)
        self.log.debug("Finished Integration")
        return ys, ms, data
        
    def surface(self):
        """Run an integration from the surface to the inner edge"""
        import scipy.integrate
        self.log.debug("Getting Outer Boundaries")
        outer_ic = outer_boundary(R=self.R_Guess,L=self.L_Guess,M=self.M,mu=self.mu,optable=self.opacity)
        ms = np.logspace(31,np.log10(self.M-1e28),5)[::-1]
        self.log.debug("Starting Integration")
        ys, data = scipy.integrate.odeint(self.integral,outer_ic,ms,full_output=True)
        self.log.debug("Finished Integration")
        return ys, ms, data
    
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
        