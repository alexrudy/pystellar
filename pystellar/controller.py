# -*- coding: utf-8 -*-
# 
#  controller.py
#  pystellar
#  
#  Created by Jaberwocky on 2012-11-09.
#  Copyright 2012 Jaberwocky. All rights reserved.
# 

from .threading import ObjectThread
from .star import Star
from .dashboard import Dashboard
from .newton import NRSolver

import logging, logging.config, logging.handlers
import time

import numpy as np

from pyshell.base import CLIEngine
from textwrap import fill


class StarEngine(CLIEngine):
    """A stellar model engine"""
    
    defaultcfg = "Star.yml"
    
    module = __name__
    
    @property
    def description(self):
        """Text description of the star engine."""
        return fill(u"Stellar Modelling in Python.")
    
    def __init__(self):
        super(StarEngine, self).__init__()
        self._threads = {}
        self._stars = {}
        self.log = logging.getLogger(__name__)
        self.datalog = logging.getLogger('data')
        self.datalog.propagate = False
        self._parser.add_argument('-q',action='store_true',
            dest='quiet', help="Run quietly")
        self._parser.add_argument('-v',action='store_true',
            dest='verbose', help="Enable verbosity.")
        
    def configure(self):
        """Configure the engine"""
        super(StarEngine, self).configure()
        if self.opts.quiet:
            self.config["Logging.handlers.console.level"] = logging.WARNING
        if self.opts.verbose:
            self.config["Logging.handlers.console.level"] = logging.DEBUG
        logging.captureWarnings(True)
        logging.config.dictConfig(self.config["Logging"])
        
    def parse(self):
        """Parse arguments which directly control the system."""
        self._parser.add_argument('--single', action='store_true',
            dest='single', help="Perform only a single integration")
        self._parser.add_argument('--outer', action='store_true',
            dest='single_outer', help="Perform only an outer single integration")
        self._parser.add_argument('--inner', action='store_true',
            dest='single_inner', help="Perform only an inner single integration")
        self._parser.add_argument('--inner-ic', action='store_true',
            dest='inner_ic', help="Perform only an inner initial conditions")
        self._parser.add_argument('--outer-ic', action='store_true',
            dest='outer_ic', help="Perform only an inner initial conditions")
        self._parser.add_argument('--no-scipy', action='store_const', 
            dest='integrator', const='pystellar', default='scipy', help="Use the custom integrator, not the scipy integrator.")
        self._parser.add_argument('--linear', action='store_false', 
            dest='logmode', help="Disable the logarithmic mass variable")
        self._parser.add_argument('--jac', action='store_true', 
            dest='jacobian', help="Do the jacobian")
        
        super(StarEngine, self).parse()
        self.opts.single_inner |= self.opts.single
        self.opts.single_outer |= self.opts.single
        self.opts.single = self.opts.single or self.opts.single_inner or self.opts.single_outer or self.opts.inner_ic or self.opts.outer_ic
        self.config.setdefault("System.Integrator.Method",self.opts.integrator)
        if self.opts.single:
            self.config["Dashboard.Figures.live.active"] = True
        
    def start(self):
        """Start the engine!"""
        self._start_time = time.clock()
        self._threads["master"] = Star(config=self.config.store)
        self.stars["master"] = [self.threads["master"]]
        self._threads["opacity"] = self._threads["master"].opacity
        self._threads["dashboard"] = ObjectThread(Dashboard,iargs=(self.config,),timeout=self.config["System.Dashboard.Timeout"],locking=False)
        self._threads["dashboard"].start()
        self._threads["dashboard"].create_dashboard()
        self.star_threads()
        self._threads["newton"] = ObjectThread(NRSolver,ikwargs=dict(stars=self.stars,config=self.config,dashboard=self.threads["dashboard"]),timeout=self.config["System.NewtonRapson.Timeout"],locking=True)
        self._threads["newton"].start()
        self.log.info("Threads Launched")
        if self.opts.single:
            self.run_single()
        elif self.opts.jacobian:
            self.threads["newton"].nrsolve()
            self.threads["newton"].release()
        
    def star_threads(self):
        """Launch the star threads"""
        optable_args = self._threads["opacity"].duplicator
        dashboard_args = self._threads["dashboard"].duplicator
        for star in self.config["System.Stars"]:
            self.stars[star] = []
            for n in range(self.config["System.Stars"][star]):
                star_thread_name = "%s-%d" % (star,n)
                star_thread = ObjectThread(Star,
                    ikwargs=dict(config=self.config.store,optable_args=optable_args,dashboard_args=dashboard_args),
                    timeout=self.config["System.Threading.Timeout"],locking=True,name=star_thread_name)
                self.stars[star] += [star_thread]
                self.threads[star_thread_name] = star_thread 
                star_thread.start()
    
    
    def run_single(self):
        """Operate a Single Integrator"""
        self.log.info("Starting Single Integration")
        if self.opts.single_inner:
            self.log.info("Calling Center-Integration")
            self.stars["Center"][0].center()
        if self.opts.single_outer:
            self.log.info("Calling Outer-Integration")
            self.stars["Surface"][0].surface()
        
        if self.opts.inner_ic:
            self.log.info("Calling Center-IC-Tests")
            self.stars["Center"][1].show_center_start()
            self.stars["Center"][1].release()
        if self.opts.outer_ic:
            self.log.info("Calling Surface-IC-Tests")
            self.stars["Surface"][1].show_surface_start()
            self.stars["Surface"][1].release()
            
        self.log.info("Retrieving Integration")
        if self.opts.single_inner:
            ys, ms, data  = self.stars["Center"][0].retrieve()
            self.log.info("Retrieved Inner Integration")
        if self.opts.single_outer:
            ys, ms, data  = self.stars["Surface"][0].retrieve()
            self.log.info("Retrieved Outer Integration")
        self.dashboard.save("live","Single_Integration.pdf")
            
    def end(self):
        """Things to do at the end of every run!"""

        total_time = time.clock()-self._start_time
        self.log.info("Total time taken: %fs" % total_time)
        for thread in self.threads:
            if thread is not "dashboard":
                self.threads[thread].stop()
        raw_input("To End the Program, press [enter]...\n")
        self.log.info("Ending Dashboard Process")
        self.dashboard.stop()
        
    @property
    def threads(self):
        """Access for threads"""
        return self._threads
        
    @property
    def stars(self):
        """Access for threads"""
        return self._stars
        
    @property
    def dashboard(self):
        """docstring for dashboard"""
        return self.threads["dashboard"]
        
    def kill(self):
        """docstring for kill"""
        for star in self._threads.values():
            star.kill()
    
if __name__ == '__main__':
    print("Running from file: {arg}".format(arg=sys.argv[0]))
    ENGINE = StarEngine()
    ENGINE.arguments()
    ENGINE.run()