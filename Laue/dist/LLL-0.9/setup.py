#!/usr/bin/env python
from distutils.core import setup

setup(name = "LLL",
      version = "0.9",
      description = "LLL - libraries for Laue lens calculation",
      author = "Alessandro Pisa",
      author_email = "alessandro...pisa@@@gmail...com",
      url = "http://darkmoon.altervista.org",
      packages = ["LLL"],
      package_dir = {'LLL': 'LLL',},
      package_data = {'LLL': ['data/*.dat'],},
      long_description = """
LLL - libraries for Laue lens calculation
""",
)
