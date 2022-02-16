.. _directories:

Important Directories
========================

On GitHub
''''''''''''''''''''''''''''''

* **PACMAN/run**

 - **PACMAN/run/run_pacman.py**

This python script runs PACMAN. As the installation with setup.py is not tested yet, the user has to modify ``line 2`` inside of this file before being able to run the code:

.. code-block:: python

	sys.path.append('/home/zieba/Desktop/Projects/Open_source/PACMAN/')

The path in ``sys.path.append`` has to be changed to the location of PACMAN on the user's device.

    .. note:: The path should point to ``/PACMAN/`` and not to ``/PACMAN/pacman/``.

The user can now run any stage by adding it as an argument in the terminal. E.g., for running Stage 00:

.. code-block:: console

	python run_pacman.py --s00

To get more information, type:

.. code-block:: console

	python run_pacman.py --help


 - **PACMAN/run/fit_par.txt**

The fit_par file is used in Stage 30 to fit the light curve. The user defines in here which fit parameters should be fixed, shared accross visits, and sets other information like priors.


 - **PACMAN/run/obs_par.pcf**

The PACMAN control file (pcf): the user sets here which plots should be saved, the path to the data and many other parameters. A thorough explanation of all the parameters in the pcf can be found on Read The Docs: :ref:`pcf`.



* **PACMAN/pacman**

All code PACMAN needs to run, is stored here.


 - **PACMAN/pacman/reduction**

This directory contains the main scripts for the individual stages. The scripts here use many files which are saved in ``PACMAN/pacman/lib``.


 - **PACMAN/pacman/lib**

This directory contains auxiliary scripts for the stages. E.g., ``PACMAN/pacman/lib/plots.py`` creates and saves plots.


 - **PACMAN/pacman/ancil**


  + **PACMAN/pacman/ancil/bandpass**

This directory contains the bandpass of the G102 and the G141 grisms. These files will be used in Stage 03 to create the reference spectrum.


  + **PACMAN/pacman/ancil/flats**

The flats (for G102 and G141) are used to find the locations of bad pixels.


  + **PACMAN/pacman/ancil/stellar_models**

This directory contains information for PACMAN which stellar models are available to download.
PACMAN offers the user to download three different stellar models from the internet: Kurucz stellar models 1993, Castelli and Kurucz stellar models 2004 and Phoenix models by Allard and collaborators.
These models will be used in Stage 03 to create the reference spectrum.
More on this `further down <https://pacmandocs.readthedocs.io/en/latest/quickstart.html#stage-03>`_ at the walkthrough of Stage 03.


When running PACMAN
''''''''''''''''''''''''''''''

* **run directory**:

Contents:

 - run_pacman.py

 - obs_par.pcf

 - fit_par

Example: ``/home/zieba/Desktop/Projects/Open_source/PACMAN/run``.

    .. note:: | The pcf file in the run directory is ONLY used in Stage 00. It will be copied over to the work directory. The copied pcf file in the work directory will then be the pcf file for all following stages. The same is true for the fit_par.txt file. So, after running Stage 00, PACMAN does not care anymore about the changes made to the pcf file and the fit_par file in the run directory.


* **work directory**:

This directory will be created when running Stage 00.
All the results of the following stages will be stored here.

Example: ``/home/zieba/Desktop/Projects/Open_source/PACMAN/run/run_2022-01-19_16-46-19_GJ1214_Hubble13021``.
It therefore has the following form:

.. code-block:: python

    datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
    meta.workdir = 'run_' + datetime + '_' + meta.eventlabel


* **data directory**:

Example: ``/home/zieba/Desktop/Data/GJ1214_Hubble13021``.

This directory should contain the .fits files which will be reduced and analyzed.


* **pipeline directory**:

This is the heart of PACMAN containing all the code to run the different Stages.

Example: ``/home/zieba/Desktop/Projects/Open_source/PACMAN/pacman``
