#!/usr/bin/env python
import os
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    REQUIRES = f.read().splitlines()

FILES = []
for root, _, files in os.walk("PACMAN"):
    FILES += [os.path.join(root.replace("PACMAN/", ""), fname) \
        for fname in files if not fname.endswith(".py") and not fname.endswith(".pyc")]

print(FILES)

setup(name='PACMAN',
      version='0.0.1',
      description='Pipeline for reduction and analysis of HST/WFC3 G102 and G141 data',
      packages=find_packages(".", exclude=["*.tests"]),
      package_data={'PACMAN': FILES},
      install_requires=REQUIRES,
      author='Sebastian Zieba',
      author_email='zieba@mpia.de',
      url='https://github.com/sebastian-zieba/PACMAN',
      license='MIT',
      long_description='',
      zip_safe=True,
      use_2to3=False
      )
