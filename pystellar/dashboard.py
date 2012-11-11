# -*- coding: utf-8 -*-
# 
#  dashboard.py
#  pystellar
#  
#  Created by Jaberwocky on 2012-11-09.
#  Copyright 2012 Jaberwocky. All rights reserved.
# 

import numpy as np

class Dashboard(object):
    """An object for managing a pyplot dashboard."""
    def __init__(self):
        super(Dashboard, self).__init__()
        import matplotlib.pyplot as plt
        self.plt = plt
        self.plt.ion()
        self.figures = {}
        self.axes = {}
        self.lines = {}
        self.keys = ["radius","luminosity","pressure","temperature","density"]
        for key in self.keys:
            self.lines[key] = {}
        
    def create_dashboard(self):
        """Create the correct items for dashboards."""
        self.figures["main"] = self.plt.figure(figsize=(15,10))
        # self.figures["main"].canvas.resize(1000,1000)
        self.axes["radius"] = self.figures["main"].add_subplot(3,2,1)
        self.axes["radius"].set_xlabel("Mass (g)")
        self.axes["radius"].set_ylabel("Radius (cm)")
        self.axes["luminosity"] = self.figures["main"].add_subplot(3,2,2)
        self.axes["luminosity"].set_xlabel("Mass (g)")
        self.axes["luminosity"].set_ylabel("Luminosity")
        self.axes["pressure"] = self.figures["main"].add_subplot(3,2,3)
        self.axes["pressure"].set_xlabel("Mass (g)")
        self.axes["pressure"].set_ylabel(r"Pressure (Dyne/$\textrm{cm}^2$)")
        self.axes["temperature"] = self.figures["main"].add_subplot(3,2,4)
        self.axes["temperature"].set_xlabel("Mass (g)")
        self.axes["temperature"].set_ylabel("Temperature (K)")
        self.axes["density"] = self.figures["main"].add_subplot(3,2,5)
        self.axes["density"].set_xlabel("Mass (g)")
        self.axes["density"].set_ylabel(r"Density ($\textrm{g}/\textrm{cm}^3$)")
        
    def insert_data(self,xs,ys,ln):
        """docstring for insert_data"""
        if ln not in self.lines["radius"]:
             self.add_lines(xs,ys,ln)
        else:
             self.update_lines(xs,ys,ln)
        
    def append_data(self,xs,ys,ln):
        """docstring for append_data"""
        self.add_lines(xs,ys,ln)
        
    def add_lines(self,xs,ys,ln):
        """Add data to the dashboard plots"""
        m = np.asarray(xs)
        r,l,P,T = np.asarray(ys).T
        self.lines["radius"][ln], = self.axes["radius"].semilogx(m,r,"bo-")
        self.lines["luminosity"][ln], = self.axes["luminosity"].semilogx(m,l,"go-")
        self.lines["pressure"][ln], = self.axes["pressure"].semilogx(m,P,"ro-")
        self.lines["temperature"][ln], = self.axes["temperature"].semilogx(m,T,"co-")
        
    def add_density(self,xs,ys,ln):
        """docstring for add_density"""
        m = np.asarray(xs)
        rho = np.asarray(ys)
        self.lines["density"][ln], = self.axes["density"].semilogx(m,rho,"bo-")
        
    def update_lines(self,xs,ys,ln):
        """Update the dashboard plots"""
        m = np.asarray(xs)
        r,l,P,T = np.asarray(ys).T
        self.lines["radius"][ln].set_ydata(r)
        self.lines["luminosity"][ln].set_ydata(l)
        self.lines["pressure"][ln].set_ydata(P)
        self.lines["temperature"][ln].set_ydata(T)
        
    def update(self):
        """Update the drawing"""
        self.figures["main"].canvas.draw()
        
        