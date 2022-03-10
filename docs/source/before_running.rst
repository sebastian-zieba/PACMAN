.. _before_running:

Before Running
================

Before we can run we have set up some stuff:

1) **Decide on a "run directory"**

All runs will be saved in a run directory. For the analysis of the GJ1214 data I have the following directory setup:

::

	user
	└── Desktop
		└── Projects
			└── Observations
				└── Hubble
					└── GJ1214_13021
						├── pacman_script.py
						├── fit_par.txt
						└── obs_par.pcf

In this case `GJ1214_13021` is my run directory

2) **Copy the contents of run_files into your run directory**

The contents of run_files are three different files:

 - **run_pacman.py**

 - **fit_par.txt**

 - **obs_par.pcf**

These three files are part of the package but can also be downloaded under this link: `Download here <https://downgit.github.io/#/home?url=https://github.com/sebastian-zieba/PACMAN/tree/master/pacman/run_files>`_.

3) **Set up the pcf**

Before we run Stage 00, the pcf has to be updated with information.
See the :ref:`pcf site <pcf>` for more information.
Most importantly, the user has to set the datadir and rundir in the pcf.

4) Run Stage 00

You run PACMAN by using the pacman_scipt.py file. For more information see :ref:`Stage 00 <stage00>`.
