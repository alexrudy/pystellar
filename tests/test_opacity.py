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
        """Interpolated Values Succeed. logrho=[0.3,-4.0], logT=[6.3,5.0]"""
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
        """Composition should fail if X+Y > 1.0."""
        self.o.composition(X=0.7,Y=0.7)
        
    @nt.raises(ValueError)
    def test_lower_logR_bound(self):
        """Values below logR=-8.0 fail."""
        self.o.composition(X=0.70,Y=0.28)
        logrho, logT = self.o.invert_points(logR=-9.0,logT=5.00).T
        v1 = self.o.lookup(logrho=logrho,logT=logT)
        
    @nt.raises(ValueError)
    def test_upper_logR_bound(self):
        """Values above logR=1.0 fail."""
        self.o.composition(X=0.70,Y=0.28)
        logrho, logT = self.o.invert_points(logR=2.0,logT=5.00).T
        v1 = self.o.lookup(logrho=logrho,logT=logT)
        
    @nt.raises(ValueError)
    def test_lower_logT_bound(self):
        """Values below logT=3.00 fail."""
        self.o.composition(X=0.70,Y=0.28)
        logrho, logT = self.o.invert_points(logR=-4.0,logT=3.00).T
        v1 = self.o.lookup(logrho=logrho,logT=logT)
        
    @nt.raises(ValueError)
    def test_upper_logT_bound(self):
        """Values above logT=9.00 fail."""
        self.o.composition(X=0.70,Y=0.28)
        logrho, logT = self.o.invert_points(logR=-4.0,logT=9.00).T
        v1 = self.o.lookup(logrho=logrho,logT=logT)
        
    @nt.raises(ValueError)
    def test_corner_a_NaNs(self):
        """Values in the bottom right corner of the table should fail. logR=0.5, logT=8.5"""
        self.o.composition(X=0.70,Y=0.28)
        logrho, logT = self.o.invert_points(logR=0.5,logT=8.50).T
        v1 = self.o.lookup(logrho=logrho,logT=logT)
        assert v1 == 1.0
        
    @nt.raises(ValueError)
    def test_corner_b_NaNs(self):
        """Values in the top left corner of the Y=1.000 table should fail. logR=-8.0, logT=3.75"""
        self.o.composition(X=0.00,Y=1.00)
        logrho, logT = self.o.invert_points(logR=-8.0,logT=3.75).T
        v1 = self.o.lookup(logrho=logrho,logT=logT)
        assert v1 == 1.0
        
    def test_corner_a_valid(self):
        """Values in the bottom left corner of the table should succeed. logR=-8.0, logT=8.70"""
        self.o.composition(X=0.70,Y=0.28)
        logrho, logT = self.o.invert_points(logR=-8.0,logT=8.70).T
        v1 = self.o.lookup(logrho=logrho,logT=logT)
        a1 = -0.582
        assert v1 - a1  < 1e-10, u"κ mismatch: %g ≠ %g" % (v1, a1)
        
    def test_corner_b_valid(self):
        """Values in the top right corner of the table should succeed. logR=1.0, logT=3.75"""
        self.o.composition(X=0.70,Y=0.28)
        logrho, logT = self.o.invert_points(logR=1.0,logT=3.75).T
        v1 = self.o.lookup(logrho=logrho,logT=logT)
        a1 = 0.131
        assert v1 - a1  < 1e-10, u"κ mismatch: %g ≠ %g" % (v1, a1)
        
    def test_midtable_valid(self):
        """Direct values succeed. logR=-4.0, logT=5.45"""
        self.o.composition(X=0.70,Y=0.28)
        logrho, logT = self.o.invert_points(logR=-4.0,logT=5.45).T
        v1 = self.o.lookup(logrho=logrho,logT=logT)
        a1 = 0.680
        assert v1 - a1  < 1e-10, u"κ mismatch: %g ≠ %g" % (v1, a1)
