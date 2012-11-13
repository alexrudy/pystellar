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
    def __init__(self,config):
        super(Dashboard, self).__init__()
        self.config = config
        self.mpl_setup() # Performs MPL setup commands
        import matplotlib.pyplot as plt
        self.plt = plt
        self.log = logging.getLogger(__name__)
        self.figures = {}
        self.axes = {}
        self.lines = {}
        self.active_figures = []
    
    def mpl_setup(self):
        """Change matplotlib defaults as configuration requires"""
        import matplotlib as mpl
        mpl.rcParams['text.usetex'] = self.config["System.Matplotlib.text.usetex"]
        
    def create_dashboard(self):
        """Create the correct items for dashboards."""
        for figure in self.config["Dashboard.Figures"]:
            x_sub, y_sub = self.config["Dashboard.Figures"][figure]["layout"]
            self.figures[figure] = self.plt.figure(figsize=tuple(self.config["Dashboard.Figures"][figure]["size"]))
            self.figures[figure].suptitle(self.config["Dashboard.Figures"][figure].get("title",figure))
            self.axes[figure] = {}
            self.lines[figure] = {}
            for n,axes in enumerate(self.config["Dashboard.Figures"][figure]["axes"]):
                nn = self.config["Dashboard.Figures"][figure]["axes"][axes].get("n",n+1)
                self.axes[figure][axes] = self.figures[figure].add_subplot(x_sub,y_sub,nn)
                self.axes[figure][axes].set_xlabel(self.config["Dashboard.Figures"][figure]["axes"][axes]["x"])
                self.axes[figure][axes].set_ylabel(self.config["Dashboard.Figures"][figure]["axes"][axes]["y"])
                self.lines[figure][axes] = {}
            if self.config["Dashboard.Figures"][figure].get("active",True):
                self.active_figures += [figure]
        
    def add_data(self,x,y,figure,axes,line,method='plot',**kwargs):
        """Add data explicitly to a figure and axes, with name line"""
        lines = getattr(self.axes[figure][axes],method)(x,y,**kwargs)
        n = 0
        for lineart in lines:
            if line in self.lines[figure][axes]:
                self.lines[figure][axes][line].remove()
                del self.lines[figure][axes][line]
            self.lines[figure][axes][line] = lineart
            n += 1
            line += "%2d" % n
        return self.lines[figure][axes].keys()
        
    def update_data(self,x,y,figure,axes,line,**kwargs):
        """Update an existing line's data to be the new data passed in."""
        if line not in self.lines[figure][axes]:
            self.add_data(x,y,figure,axes,line,**kwargs)
        self.lines[figure][axes][line].set_data(x,y)
        
    def append_data(self,x,y,figure,axes,line,**kwargs):
        """Append data to an existing line."""
        if line not in self.lines[figure][axes]:
            self.add_data(x,y,figure,axes,line,**kwargs)
        data_x,data_y = self.lines[figure][axes][line].get_data()
        new_x = np.hstack((data_x,x))
        new_y = np.hstack((data_y,y))  
        self.lines[figure][axes][line].set_data(new_x,new_y)
        
    
    def update(self,*figures):
        """Update the drawing"""
        if len(figures) == 0:
            figures = tuple(self.figures.keys())
        for figure in figures:
            for ax in self.axes[figure].values():
                ax.relim()
                ax.autoscale_view()
            if figure in self.active_figures:
                self.figures[figure].show()
                self.figures[figure].canvas.draw()
        self.log.debug("Dashboard Updated")
        
    def save(self,figure,filename):
        """Save the final integration result"""
        self.figures[figure].savefig(filename)
        
        