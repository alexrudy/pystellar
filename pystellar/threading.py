# -*- coding: utf-8 -*-
# 
#  threading.py
#  pystellar
#  
#  Created by Jaberwocky on 2012-10-12.
#  Copyright 2012 Jaberwocky. All rights reserved.
# 

from __future__ import division

from StringIO import StringIO

import numpy as np
import scipy as sp

from pkg_resources import resource_filename
from warnings import warn

from multiprocessing import Queue, Process, Pool, cpu_count

import logging
# Alex's modules
from AstroObject.config import DottedConfiguration

class ThreadStateError(Exception):
    """Handles errors with this module's state."""
    pass
        

class ObjectsThread(object):
    """A manager for handling many threaded objects which take from a single job queue and return to a single job queue."""
    def __init__(self, Oclass, *args, **kwargs):
        super(ObjectsThread, self).__init__()
        self.Oclass = Oclass
        self._args = args
        self._nprocs = kwargs.pop('nprocs',cpu_count())
        self._kwargs = kwargs
        self._started = False
        self.input = Queue()
        self.output = Queue()
        self._kwargs['input_Q'] = self.input
        self._kwargs['output_Q'] = self.output
        self._procs = [ObjectThread(self.Oclass,*self._args,**self._kwargs) for i in xrange(self._nprocs)]
        
    @property
    def started(self):
        return self._started
        
    def start(self):
        """docstring for fname"""
        if self.started:
            raise ThreadStateError("Thread pool already started")
        for proc in self._procs:
            proc.start()
        self._started = True
        
    def __getattr__(self,attr):
        """Call a method on the underlying threaded object"""
        
        def method(*args,**kwargs):
            """A threaded method"""
            self.input.put((attr,args,kwargs))
            
        return method
    
    @property
    def empty(self):
        """Lets us know if the output queue is empty"""
        return self.output.empty()
    
    def retrieve(self,inputs=False):
        """Retrieve a return value off the top of the output queue"""
        func,args,kwargs,rvalue = self.output.get(timeout=10)
        if inputs:
            return func,args,kwargs,rvalue
        else:
            return rvalue
        
    def stop(self):
        """Send the thread stop signal."""
        for proc in self._procs:
            self.input.put((proc.STOP,None,None),timeout=10)
        while not self.output.empty():
            self.retrieve()
        for proc in self._procs:
            proc.join(timeout=10)

        
    

class ObjectThread(Process):
    """A thread management object for single object threads"""
    def __init__(self, Oclass, *args, **kwargs):
        super(ObjectThread, self).__init__()
        self.Oclass = Oclass
        self.input = kwargs.pop('input_Q',Queue())
        self.output = kwargs.pop('output_Q',Queue())
        self._args = args
        self._kwargs = kwargs
        
    STOP = '..stop'
    
    def __getattr__(self,attr):
        """Call a method on the underlying threaded object"""
        
        def method(*args,**kwargs):
            """A threaded method"""
            self.input.put((attr,args,kwargs),timeout=10)
            
        return method
    
    def retrieve(self,inputs=False):
        """Retrieve a return value off the top of the output queue"""
        func,args,kwargs,rvalue = self.output.get(timeout=10)
        if inputs:
            return func,args,kwargs,rvalue
        else:
            return rvalue
        
    def stop(self):
        """Send the thread stop signal."""
        self.input.put((self.STOP,None,None),timeout=10)
        while not self.output.empty():
            self.retrieve()
        self.join()
        
    def run(self):
        """starts an opactiy thread which takes a queue of items."""
        O = self.Oclass(*self._args, **self._kwargs)
        done = False
        while not done:
            func, args, kwargs = self.input.get(timeout=10)
            if self.STOP == func:
                done = True
            else:
                try:
                    attr = getattr(O,func)
                    if callable(attr):
                        rvalue = attr(*args,**kwargs)
                    elif len(args)==0 and len(kwargs)==0:
                        rvalue = attr
                    else:
                        raise AttributeError("Asked for attribute with arguments!")
                    if rvalue is not None:
                        self.output.put((func,args,kwargs,rvalue),timeout=10)
                except Exception:
                    done = True
                    raise
