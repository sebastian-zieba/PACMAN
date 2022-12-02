.. _directories:

Repository Structure
========================

In the package
''''''''''''''''''''''''''''''


* `PACMAN/src/pacman <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman>`_

  All code PACMAN needs to run and templates for the pcf and fit_par files are stored here.
  The scripts here use other python code from ``PACMAN/src/pacman/lib``.

  - `PACMAN/src/pacman/lib <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/lib>`_

  This directory contains auxiliary scripts for the stages.
  E.g., ``PACMAN/src/pacman/lib/plots.py`` creates and saves plots.


  - `PACMAN/src/pacman/data <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data>`_

    + `PACMAN/src/pacman/data/bandpass <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/bandpass>`_

    This directory contains the throughput for the G102 and the G141 grisms.
    These files will be used in Stage 03 to create the reference spectrum.


    + `PACMAN/src/pacman/data/flats <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/flats>`_

    The flats (for the G102 and G141 grisms) are used to find the locations of bad pixels.


    + `PACMAN/src/pacman/data/stellar_models <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/stellar_models>`_

    This directory tells PACMAN which stellar models are available to download.
    PACMAN can download three different stellar models: `Kurucz stellar models 1993 <https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/kurucz-1993-models>`_, `Castelli and Kurucz stellar models 2004 <https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas>`_ and `Phoenix models by Allard and collaborators <https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/phoenix-models-available-in-synphot>`_.
    These models will be used in Stage 03 to create the reference spectrum.
    More on this at the walkthrough of `Stage 03 <https://pacmandocs.readthedocs.io/en/latest/quickstart.html#stage-03>`_ using GJ1214 data as an example.


    + `PACMAN/src/pacman/data/run_files <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files>`_

    .. note:: These three files have to be copied over to the run directory where the analysis will be run from! You find more information about that in :ref:`Before Running <before_running>`.


    * -> `PACMAN/src/pacman/data/run_files/pacman_script.py <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files/pacman_script.py>`_

    This python script runs PACMAN. If the user did not install PACMAN properly as explained in the :ref:`installation walkthrough <installation>` (not recommended) he or she might have to add a path to the sys.path.append pointing to the PACMAN directory:

    .. code-block:: python

  	  sys.path.append('/home/zieba/Desktop/Projects/Open_source/PACMAN/')

    The path in ``sys.path.append`` has to be changed to the location of PACMAN on the user's device.

    .. note:: The path should point to ``/PACMAN/`` and not to ``/PACMAN/src/pacman/``.


    * -> `PACMAN/src/pacman/data/run_files/obs_par.pcf <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files/obs_par.pcf>`_

    The PACMAN control file (pcf): the user sets parameters related to the analysis. E.g., which plots should be saved, the path to the data directory and others.
    A thorough explanation of all the parameters in the pcf can be found on the :ref:`PCF page<pcf>`.


    * -> `PACMAN/src/pacman/data/run_files/fit_par.txt <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files/fit_par.txt>`_

    The fit_par.txt file is used in Stage 30 to fit the light curve. The user defines here which fit parameters should be fixed, shared across visits, and sets other information like priors.




Nomenclature
''''''''''''''''''''''''''''''

An example for a directory structure:

	user
	└── Desktop
		└── Projects
			└── Observations
				└── Hubble
					└── "run directory"
						├── pacman_script.py
						├── fit_par.txt
						├── obs_par.pcf
						└── "work directory"

		└── "data directory"
			└── GJ1214_Hubble13021
                            ├── ibxy06d0q_ima.fits
                            ├── ...
                            └── ibxy07ryq_ima.fits



* **run directory**:

  Example: ``/home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021``.

  .. note:: You have to copy the contents of run_files into this directory.

  Contents to copy into the run directory:

   - pacman_script.py

   - obs_par.pcf

   - fit_par.txt

   If you pip-installed, downloaded or cloned the GitHub repository, you'll find the run_files directory (with templates for these three files) in ``PACMAN/src/pacman/data/run_files``.
   They can also be downloaded under this link: `Download here <https://downgit.github.io/#/home?url=https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files>`_.
   You have to copy these files into your run directory.

    .. note:: | The pcf file in the run directory is ONLY used in Stage 00.
              | When running Stage 00, the pcf and fit_par files will be copied over to the work directory.
              | The copied pcf file in the work directory will then be the pcf file for all following stages.
              | The same is true for the fit_par.txt file.
              | So, after running Stage 00, PACMAN does not care anymore about the changes made to the pcf file and the fit_par file in the run directory!


* **work directory**:

  Example: ``/home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021/run_2022-03-04_15-10-29_GJ1214_Hubble13021``.

  This directory will be created in Stage 00.
  All the results of the following stages will be stored here.

  The name of the work directory is a combination of the following parts:

  "run_" + "YYYY-MM-DD_HH-MM-SS_" + "eventlabel"

  So for example: run_2022-03-04_15-10-29_GJ1214_Hubble13021

  The eventlabel is chosen by the user when running Stage 00.


* **data directory**:

  Example: ``/home/zieba/Desktop/Data/GJ1214_Hubble13021``.

  This directory should contain the .fits files which will be reduced and analyzed.


* **pipeline directory**:

  Example: ``/home/zieba/Desktop/Projects/Open_source/PACMAN/src``

  This is the heart of PACMAN containing all the code and data to run the different stages.
