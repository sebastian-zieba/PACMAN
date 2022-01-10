# #!/usr/bin/env python
# import os
# from setuptools import setup, find_packages
#
# REQUIRES = ['batman-package',
#             'h5py',
#             'ipython',
#             'matplotlib',
#             'numpy>=1.20.0',
#             'pytest',
#             'pyyaml',
#             'scipy']
#
# DEPENDENCY_LINKS = []
#
# FILES = []
# for root, _, files in os.walk("wfc3-pipeline"):
#     FILES += [os.path.join(root.replace("wfc3-pipeline/", ""), fname) \
#         for fname in files if not fname.endswith(".py") and not fname.endswith(".pyc")]
#
# setup(name='wfc3-pipeline',
#       version='0.0.1',
#       description='Lightcurve fitting package for time-series observations',
#       packages=find_packages(".", exclude=["*.tests"]),
#       package_data={'wfc3-pipeline': FILES},
#       install_requires=REQUIRES,
#       dependency_links=DEPENDENCY_LINKS,
#       #author='Section 5',
#       #author_email='kbstevenson@gmail.com',
#       license='MIT',
#       #url='https://github.com/kevin218/Eureka',
#       long_description='',
#       zip_safe=True,
#       use_2to3=False
# )


from setuptools import setup, find_packages
setup(name='PACMAN',
      version='0.0.1',
      packages=find_packages(".")
      )
