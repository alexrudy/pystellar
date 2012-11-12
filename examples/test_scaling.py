#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  test_scaling.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-11-11.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 


import numpy as np
import matplotlib.pyplot as plt

M = 33
fp = 28
x = np.linspace(28,33,50)
m = np.power(10,M - x)

print m