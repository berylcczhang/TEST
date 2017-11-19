#!/usr/bin/env python

import os
import os.path as osp
import shutil
import sys
from distutils.core import setup, Extension
import subprocess
import numpy
if 'develop' in sys.argv:
    # use setuptools for develop, but nothing else
    from setuptools import setup

with open('README.md') as file:
    long_description = file.read()

with open('CHANGES') as file:
    long_description += file.read()

with open('REQUIREMENTS') as file:
    requirements = file.readlines()

try:  # Python 3.x
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:  # Python 2.x
    from distutils.command.build_py import build_py

from slugpy import __version__ as version
tagname = "slugpy_%s" % (version)

download_url = "https://bitbucket.org/krumholz/slug2/"

bpc = osp.join('slugpy', 'bayesphot', 'bayesphot_c')

setup(name='SLUGPY',
      version=version,
      description='a python package to use with the SLUG stochastic SPS code',
      long_description=long_description,
      author=['Mark Krumholz'],
      author_email=['mark.krumholz@gmail.com'],
      packages=['slugpy', 'slugpy.bayesphot', 'slugpy.cloudy', 
                'slugpy.cluster_slug', 'slugpy.sfr_slug'],
      requires=requirements,
      cmdclass={'build_py': build_py},
      classifiers=[
                   "Development Status :: 3 - Alpha",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
                  ],
      ext_modules = [Extension('slugpy.bayesphot.bayesphot',
                               [osp.join(bpc, 'geometry.c'), 
                                osp.join(bpc, 'kdtree.c'),
                                osp.join(bpc, 'kernel_density.c'),
                                osp.join(bpc, 'kernel_density_neighbors.c'),
                                osp.join(bpc, 'kernel_density_rep.c'),
                                osp.join(bpc, 'kernel_density_util.c'),
                                osp.join(bpc, 'kernel_density_draw.c')
                            ],
                               libraries=['gsl', 'gslcblas'])]
  )
