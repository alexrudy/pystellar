# -*- coding: utf-8 -*-
# 
#  controller.py
#  pystellar
#  
#  Created by Jaberwocky on 2012-11-09.
#  Copyright 2012 Jaberwocky. All rights reserved.
# 

from pystellar.threading import ObjectThread
from pystellar.star import Star

import logging

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
        self._stars = {}
        self.log = logging.getLogger(__name__)
        self.__setup_main_log()
        
    def __setup_main_log(self):
        """Setup the main pystellar logger"""
        self._main_log = logging.getLogger("pystellar")
        self._stream = logging.StreamHandler()
        self._stream.setFormatter(logging.Formatter(fmt="--> %(message)s"))
        self._main_log.addHandler(self._stream)
        self._main_log.setLevel(logging.DEBUG)
        
    def parse(self):
        """docstring for parse"""
        self._parser.add_argument('--single', action='store_true',
            dest='single', help="Perform only a single integration")
        super(StarEngine, self).parse()
        
    def start(self):
        """Start the engine!"""
        self._stars["master"] = Star(filename=self._opts.config)
        if self._opts.single:
            self.run_single()
        
        
    
    def run_single(self):
        """docstring for run_single"""
        self.log.info("Starting Single Integration")
        optable_args = self._stars["master"].opacity.duplicator
        self._stars["inner_1"] = ObjectThread(Star,ikwargs=dict(filename=self._opts.config,optable_args=optable_args),timeout=120)
        self._stars["inner_1"].start()
        self._stars["outer_1"] = ObjectThread(Star,ikwargs=dict(filename=self._opts.config,optable_args=optable_args),timeout=120)
        self._stars["outer_1"].start()
        
        self.log.info("Calling Center-Integration")
        self._stars["inner_1"].center()
        self.log.info("Calling Outer-Integration")
        self._stars["outer_1"].surface()
        
    def end(self):
        """docstring for end"""
        if self._opts.single:
            self.log.info("Retrieving Single Integration")
            ys, ms, data  = self._stars["inner_1"].retrieve()
            import matplotlib.pyplot as plt
            # plt.ion()
            fig = plt.figure()
            ax1 = fig.add_subplot(221)
            ax1.loglog(ms,ys[:,0],".-")
            ax1.set_xlabel("Mass")
            ax1.set_ylabel("Radius")
            ax2 = fig.add_subplot(222)
            ax2.loglog(ms,ys[:,1],".-")
            ax2.set_xlabel("Mass")
            ax2.set_ylabel("Luminosity")
            ax3 = fig.add_subplot(223)
            ax3.loglog(ms,ys[:,2],".-")
            ax3.set_xlabel("Mass")
            ax3.set_ylabel("Pressure")
            ax4 = fig.add_subplot(224)
            ax4.loglog(ms,ys[:,3],".-")
            ax4.set_xlabel("Mass")
            ax4.set_ylabel("Temperature")
            plt.show()
            ys, ms, data  = self._stars["outer_1"].retrieve()
            ax1.loglog(ms,ys[:,0],".-")
            ax2.loglog(ms,ys[:,1],".-")
            ax3.loglog(ms,ys[:,2],".-")
            ax4.loglog(ms,ys[:,3],".-")
        for star in self._stars.values():
            star.stop()
        
    def kill(self):
        """docstring for kill"""
        for star in self._stars.values():
            star.kill()
    
if __name__ == '__main__':
    print("Running from file: {arg}".format(arg=sys.argv[0]))
    ENGINE = StarEngine()
    ENGINE.arguments()
    ENGINE.run()