#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  thread_opacity.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-12.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

from __future__ import division

from pystellar.opacity import OpacityThread

print "--Launching Thread"
OT = OpacityThread(fkey='GN93hz')
print "--Thread Launched"
OT.composition(X=0.70,Y=0.28)
print "--Set Composition"
n,X,Y,Z =  OT.status()
print u"Fixed composition at X=%.3g,Y=%.3g,Z=%.3g" % (X,Y,Z)
print u"Using table %d" % (n + 1)

logTs   = [6.3, 5.0]
logrhos = [0.3,-4.0]

print "--Doing Lookups"
for logT,logrho in zip(logTs,logrhos):
    OT.lookup(logT=logT,logrho=logrho)

print "--Doing Retrieves"
for i in xrange(len(logTs)):
    logrho, logT, logK = OT.retrieve()
    print u"log(T)=%5.3f, log(ρ)=%6.3f\n → log(κ)=%5.3f" % (logT,logrho,logK)
    
OT.stop()
print "Threads Done!"