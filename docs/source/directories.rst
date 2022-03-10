.. _directories:

Main Directories
========================

.. image:: media/20220310_113633.jpg

In the package
''''''''''''''''''''''''''''''


* **PACMAN/pacman**

  All code PACMAN needs to run and templates for the pcf and fit_par files are stored here.


  - **PACMAN/pacman/ancil**

    + **PACMAN/pacman/ancil/bandpass**

    This directory contains the bandpass of the G102 and the G141 grisms.
    These files will be used in Stage 03 to create the reference spectrum.


    + **PACMAN/pacman/ancil/flats**

    The flats (for G102 and G141) are used to find the locations of bad pixels.


    + **PACMAN/pacman/ancil/stellar_models**

    This directory contains information for PACMAN which stellar models are available to download.
    PACMAN offers the user to download three different stellar models from the internet: Kurucz stellar models 1993, Castelli and Kurucz stellar models 2004 and Phoenix models by Allard and collaborators.
    These models will be used in Stage 03 to create the reference spectrum.
    More on this at the walkthrough of `Stage 03 <https://pacmandocs.readthedocs.io/en/latest/quickstart.html#stage-03>`_ using GJ1214 data as an example.


  - **PACMAN/pacman/lib**

  This directory contains auxiliary scripts for the stages.
  E.g., ``PACMAN/pacman/lib/plots.py`` creates and saves plots.


  - **PACMAN/pacman/reduction**

  This directory contains the main scripts for the individual stages.
  The scripts here use other python code which is saved in ``PACMAN/pacman/lib``.


  - **PACMAN/pacman/run_files**

    + **PACMAN/pacman/run_files/pacman_script.py**

    This python script runs PACMAN. If the user did not install pacman properly as explained in :ref:`the installation walkthrough <installation>` (not recommended) he or she might have to add a path to the sys.path.append pointing to the PACMAN directory:

    .. code-block:: python

  	  sys.path.append('/home/zieba/Desktop/Projects/Open_source/PACMAN/')

    The path in ``sys.path.append`` has to be changed to the location of PACMAN on the user's device.

    .. note:: The path should point to ``/PACMAN/`` and not to ``/PACMAN/pacman/``.

    The user can now run any stage by adding it as an argument in the terminal. E.g., for running Stage 00:

    .. code-block:: console

  	  python pacman_script.py --s00 --eventlabel='GJ1214_Hubble13021'

    To get more information, type:

    .. code-block:: console

  	  python pacman_script.py --help


    + **PACMAN/pacman/run_files/fit_par.txt**

    The fit_par file is used in Stage 30 to fit the light curve. The user defines in here which fit parameters should be fixed, shared across visits, and sets other information like priors.


    + **PACMAN/pacman/run_files/obs_par.pcf**

    The PACMAN control file (pcf): the user sets here which plots should be saved, the path to the data and many other parameters.
    A thorough explanation of all the parameters in the pcf can be found on Read The Docs: :ref:`pcf`.



Nomenclature
''''''''''''''''''''''''''''''

* **run directory**:

  Example: ``/home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021``.

  Contents:

   - pacman_script.py

   - obs_par.pcf

   - fit_par

   If you downloaded or cloned the GitHub repository it includes the run_files directory with templates for these three files.
   They can also be downloaded under this link: `Download here <https://downgit.github.io/#/home?url=https://github.com/sebastian-zieba/PACMAN/tree/master/pacman/run_files>`_.
   You have to copy these files into your run directory.

    .. note:: | The pcf file in the run directory is ONLY used in Stage 00.
              | It will be copied over to the work directory.
              | The copied pcf file in the work directory will then be the pcf file for all following stages.
              | The same is true for the fit_par.txt file.
              | So, after running Stage 00, PACMAN does not care anymore about the changes made to the pcf file and the fit_par file in the run directory.


* **work directory**:

  Example: ``/home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021/run_2022-03-04_15-10-29_GJ1214_Hubble13021``.

  This directory will be created when running Stage 00.
  All the results of the following stages will be stored here.

  It therefore has the following form:

  .. code-block:: python

      datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
      meta.workdir = 'run_' + datetime + '_' + meta.eventlabel


* **data directory**:

  Example: ``/home/zieba/Desktop/Data/GJ1214_Hubble13021``.

  This directory should contain the .fits files which will be reduced and analyzed.


* **pipeline directory**:

  Example: ``/home/zieba/Desktop/Projects/Open_source/PACMAN/pacman``

  This is the heart of PACMAN containing all the code to run the different Stages.












PACMAN consists out of these important parts:

* pacman
Contains the heart of pacman with all python scipts needed to reduce and analyse HST data.

* data directory
a local directory which contains all the ima fits files (has to be set in pcf before running Stage 00).

* run directory
a local directory where the run will be saved (has to be set in pcf before running Stage 00).
This directory will ultimately get an additional work directory every time Stage 00 is being run.

* work directory
a subdirectory of the run directory.
It has the following form: eg, ``run_2022-01-19_16-46-19_GJ1214_Hubble13021``.
It therefore contains the date and time Stage 00 has been run and the eventlabel.
