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

AstroObjectReq = "0.6"
AstroObjectDep = "AstroObject>=" + AstroObjectReq
AstroObjectVer = "0.6"
AstroObjectURL = "https://github.com/alexrudy/AstroObject/zipball/v%(ver)s#egg=AstroObject-%(ver)s" % { 'ver' : AstroObjectVer}

setup(
    name = 'pystellar',
    version = version,
    packages = find_packages(exclude=['tests']),
    package_data = {'pystellar':['OPAL.yml']},
    install_requires = ['distribute','PyYAML>=3.10','astropysics',AstroObjectDep],
    dependency_links = [AstroObjectURL],
    test_suite = 'tests',
    author = 'Alexander Rudy',
    author_email = 'dev@alexrudy.org',
    entry_points = {
        'console_scripts' : ["PyStar = pystellar.controller:StarEngine.script"]
    }
)