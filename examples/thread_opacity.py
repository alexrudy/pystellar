#!/usr/bin/env python
# -*- coding: utf-8 -*-
# #!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  thread_opacity.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-12.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

from __future__ import division

from pystellar.opacity import OpacityTable
from pystellar.threading import ObjectsThread

import numpy as np
import time

print "--Launching Thread"
start = time.clock()
OT = ObjectsThread(OpacityTable,nprocs=4,kwargs=dict(fkey='OP17',X=0.70,Y=0.28))
OT.start()
finish = time.clock()
print "--Threads Launched: %g" % (finish-start)
print "++Launching Object"
start = time.clock()
OO = OpacityTable(fkey='GN93hz',X=0.70,Y=0.28)
finish = time.clock()
print "++Object Launched: %g" % (finish-start)
print "--Get Composition"
start = time.clock()
OT.properties()
n,X,Y,Z,dXc,dXo = OT.retrieve()
print u"Fixed composition at X=%.3g,Y=%.3g,Z=%.3g" % (X,Y,Z)
print u"Using table %d" % (n + 1)
finish = time.clock()
print "--Got Composition: %g" % (finish-start)
print "++Get Composition"
start = time.clock()
n,X,Y,Z,dXc,dXo = OO.properties()
print u"Fixed composition at X=%.3g,Y=%.3g,Z=%.3g" % (X,Y,Z)
print u"Using table %d" % (n + 1)
finish = time.clock()
print "++Got Composition: %g" % (finish-start)

ntest = 1e4
logTs   = np.linspace(4.0,6.3,ntest)
logrhos = np.linspace(-6.0,0.3,ntest)
logkappa = np.empty((ntest,2))
ntest = int(ntest)

print "--Doing Interpolation"
start = time.clock()
for i in xrange(ntest):
    OT.lookup(logT=logTs[i],logrho=logrhos[i])
for i in xrange(ntest):
    func, args, kwargs, rvalue = OT.retrieve(inputs=True)
    logkappa[i,0] = rvalue
finish = time.clock()
print "--Done Interpolation: %g" % (finish-start)

print "++Doing Interpolation"
start = time.clock()
for i in xrange(ntest):
    logkappa[i,1] = OO.lookup(logT=logTs[i],logrho=logrhos[i])
finish = time.clock()
print "++Done Interpolation: %g" % (finish-start)

print "--Spurious Commands"
for i in xrange(100):
    OT.lookup(logT=logTs[i],logrho=logrhos[i])
OT.stop()
print "Threads Done!"