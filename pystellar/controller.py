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

import logging

import numpy as np

from pyshell.base import CLIEngine
from textwrap import fill


class StarEngine(CLIEngine):
    """A stellar model engine"""
    
    defaultcfg = "Star.yml"
    
    module = __name__
    
    _star_threads = ["inner_1","outer_1","inner_2","outer_2"]
    
    @property
    def description(self):
        """Text description of the star engine."""
        return fill(u"Stellar Modelling in Python.")
    
    def __init__(self):
        super(StarEngine, self).__init__()
        self._threads = {}
        self.log = logging.getLogger(__name__)
        self.__setup_main_log()
        
    def __setup_main_log(self):
        """Setup the main pystellar logger"""
        self._main_log = logging.getLogger("pystellar")
        self._stream = logging.StreamHandler()
        self._stream.setFormatter(logging.Formatter(fmt="--> %(message)s"))
        self._main_log.addHandler(self._stream)
        self._main_log.setLevel(logging.INFO)
        
    def parse(self):
        """docstring for parse"""
        self._parser.add_argument('--single', action='store_true',
            dest='single', help="Perform only a single integration")
        self._parser.add_argument('--outer', action='store_true',
            dest='outer_s', help="Perform only an outer single integration")
        self._parser.add_argument('--inner', action='store_true',
            dest='inner_s', help="Perform only an inner single integration")
        self._parser.add_argument('--no-scipy', action='store_false', 
            dest='scipy', help="Use the custom integrator, not the scipy integrator.")
        self._parser.add_argument('--no-logmode', action='store_false', 
            dest='logmode', help="Disable the logarithmic mass variable")
        self._parser.add_argument('-v',action='store_true',
            dest='verbose', help="Enable verbosity.")
        super(StarEngine, self).parse()
        
    def start(self):
        """Start the engine!"""
        if self._opts.verbose:
            self._main_log.setLevel(logging.DEBUG)
        self._threads["master"] = Star(filename=self._opts.config)
        self._threads["opacity"] = self._threads["master"].opacity
        self._threads["dashboard"] = ObjectThread(Dashboard,timeout=self.config["System.Dashboard.Timeout"],locking=False)
        self._threads["dashboard"].start()
        self._threads["dashboard"].create_dashboard()
        self._threads["dashboard"].update()
        self.star_threads()
        if self._opts.single or self._opts.inner_s or self._opts.outer_s:
            self.run_single()
        
    def star_threads(self):
        """Launch the star threads"""
        optable_args = self._threads["opacity"].duplicator
        dashboard_args = self._threads["dashboard"].duplicator
        for star in self._star_threads:
            self._threads[star] = ObjectThread(Star,
                ikwargs=dict(filename=self._opts.config,optable_args=optable_args,dashboard_args=dashboard_args),
                timeout=self.config["System.Threading.Timeout"],locking=True)
            self._threads[star].start()
            
            if self._opts.scipy:
                self._threads[star].use_scipy()
            else:
                self._threads[star].use_pystellar()
            self._threads[star].release()
            
            if self._opts.logmode:
                self._threads[star].use_logmode()
            else:
                self._threads[star].disable_logmode()
            self._threads[star].release()

    
    def run_single(self):
        """Operate a Single Integrator"""
        self.log.info("Starting Single Integration")
        if self._opts.single or self._opts.inner_s:
            self.log.info("Calling Center-Integration")
            self._threads["inner_1"].center()
        if self._opts.single or self._opts.outer_s:
            self.log.info("Calling Outer-Integration")
            self._threads["outer_1"].surface()
        
    def end(self):
        """docstring for end"""
        self.log.info("Retrieving Integration")
        if self._opts.single or self._opts.inner_s:
            ys, ms, data  = self._threads["inner_1"].retrieve()
            self.log.info("Retrieved Inner Integration")
        if self._opts.single or self._opts.outer_s:
            ys, ms, data  = self._threads["outer_1"].retrieve()
            print "NaNs Retrieved: %d / %d" % (np.sum(np.isnan(ys)),ys.size)
            self.log.info("Retrieved Outer Integration")
        raw_input("Continue...")
        for thread in self._threads.values():
            thread.stop()
        
    def kill(self):
        """docstring for kill"""
        for star in self._threads.values():
            star.kill()
    
if __name__ == '__main__':
    print("Running from file: {arg}".format(arg=sys.argv[0]))
    ENGINE = StarEngine()
    ENGINE.arguments()
    ENGINE.run()