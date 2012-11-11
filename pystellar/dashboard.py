# -*- coding: utf-8 -*-
# 
#  dashboard.py
#  pystellar
#  
#  Created by Jaberwocky on 2012-11-09.
#  Copyright 2012 Jaberwocky. All rights reserved.
# 

class Dashboard(object):
    """An object for managing a pyplot dashboard."""
    def __init__(self):
        super(Dashboard, self).__init__()
        import matplotlib.pyplot as plt
        self.plt = plt
        
    
        