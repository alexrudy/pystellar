# -*- coding: utf-8 -*-
# 
#  setup.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-03.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 
from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

from pystellar import version

setup(
    name = 'pystellar',
    version = version,
    packages = find_packages(exclude=['tests']),
    package_data = {'pystellar':['OPAL.yml']},
    install_requires = ['distribute','PyYAML>=3.10','astropysics'],
    test_suite = 'tests',
    author = 'Alexander Rudy',
    author_email = 'dev@alexrudy.org',
)