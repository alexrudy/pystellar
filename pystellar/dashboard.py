# -*- coding: utf-8 -*-
# 
#  dashboard.py
#  pystellar
#  
#  Created by Jaberwocky on 2012-11-09.
#  Copyright 2012 Jaberwocky. All rights reserved.
# 

import numpy as np
import logging

def add_to_line(line,x,y):
    """docstring for add_to_line"""
    data_x,data_y = line.get_data()
    new_x = np.hstack((data_x,x))
    new_y = np.hstack((data_y,y))  
    line.set_data(new_x,new_y)

class Dashboard(object):
    """An object for managing a pyplot dashboard."""
    def __init__(self):
        super(Dashboard, self).__init__()
        import matplotlib.pyplot as plt
        self.plt = plt
        self.log = logging.getLogger(__name__)
        self.plt.ion()
        self.figures = {}
        self.axes = {}
        self.lines = {}
        self.keys = ["radius","luminosity","pressure","temperature","density","epsilon"]
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
        self.axes["epsilon"] = self.figures["main"].add_subplot(3,2,6)
        self.axes["epsilon"].set_xlabel("Mass (g)")
        self.axes["epsilon"].set_ylabel(r"$\epsilon$")
        
    def insert_data(self,xs,ys,ln):
        """docstring for insert_data"""
        if ln not in self.lines["radius"]:
             self.add_lines(xs,ys,ln)
        else:
             self.update_lines(xs,ys,ln)
        
    def replace_data(self,xs,ys,ln):
        """Replace data in a given line"""
        if ln not in self.lines["radius"]:
            self.add_lines(xs,ys,ln)
        else:
            m = np.asarray(xs)
            r,l,P,T = np.asarray(ys).T
            for line, data in zip(["radius","luminosity","pressure","temperature"],[r,l,P,T]):
                self.lines[line][ln].set_data(m,data)
        
    def add_lines(self,xs,ys,ln):
        """Add data to the dashboard plots"""
        m = np.asarray(xs)
        r,l,P,T = np.asarray(ys).T
        self.lines["radius"][ln], = self.axes["radius"].plot(m,r,"b-")
        self.lines["luminosity"][ln], = self.axes["luminosity"].plot(m,l,"g-")
        self.lines["pressure"][ln], = self.axes["pressure"].plot(m,P,"r-")
        self.lines["temperature"][ln], = self.axes["temperature"].plot(m,T,"c-")
        
    def add_density(self,xs,ys,ln,append=True):
        """docstring for add_density"""
        m = np.asarray(xs)
        rho = np.asarray(ys)
        if ln not in self.lines["density"]:
            self.lines["density"][ln], = self.axes["density"].plot(m,rho,"m-")
        elif append:
            add_to_line(self.lines["density"][ln],m,rho)
        else:
            self.lines["density"][ln].set_data(m,rho)
    
    def add_epsilon(self,xs,ys,ln,append=True):
        """docstring for add_density"""
        m = np.asarray(xs)
        eps = np.asarray(ys)
        if ln not in self.lines["epsilon"]:
            self.lines["epsilon"][ln], = self.axes["epsilon"].plot(m,eps,"m-")
        elif append:
            add_to_line(self.lines["epsilon"][ln],m,eps)
        else:
            self.lines["epsilon"][ln].set_data(m,eps)
            
        
    def update_lines(self,xs,ys,ln):
        """Update the dashboard plots"""
        m = np.asarray(xs)
        r,l,P,T = np.asarray(ys).T
        for line, data in zip(["radius","luminosity","pressure","temperature"],[r,l,P,T]):
            add_to_line(self.lines[line][ln],m,data)
        

    def update(self):
        """Update the drawing"""
        for ax in self.axes.values():
            ax.relim()
            ax.autoscale_view()
        self.figures["main"].canvas.draw()
        self.log.debug("Dashboard Updated")
        
    def save(self):
        """Save the final integration result"""
        self.figures["main"].savefig("Integration.pdf")
        
        