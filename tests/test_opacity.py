# -*- coding: utf-8 -*-
# 
#  test_opacity.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-04.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

import nose.tools as nt
from pystellar.opacity import OpacityTable

class test_OpacityTable(object):
    """OpacityTable"""
    
    def setup(self):
        """Set up this test suite"""
        self.o = OpacityTable("GN93hz",load=False)

    def test_solar_composition(self):
        """Test for some known values at solar composition"""
        self.o.composition(X=0.7,Y=0.28)
        assert self.o.n == 72, u"Table Mismatch, %g ≠ %g" % (self.o.n,72)
        v1 = self.o.lookup(logrho=0.3,logT=6.3)
        a1 = 1.885
        assert v1 - a1 < 1e-10, u"κ mismatch: %g ≠ %g" % (v1, a1)
        v2 = self.o.lookup(logrho=-4.0,logT=5.0)
        a2 = 3.436
        assert v2 - a2 < 1e-10, u"κ mismatch: %g ≠ %g" % (v2, a2)
        
    @nt.raises(AssertionError)
    def test_sanity_comp(self):
        """Sanity Errors for Composition"""
        self.o.composition(X=0.7,Y=0.70)
        