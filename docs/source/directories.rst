.. _directories:

Repository Structure
====================

In the package
''''''''''''''


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


    + `PACMAN/src/pacman/data/pacman_run_files <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/pacman_run_files>`_

    .. note:: These three files have to be copied over to the run directory where the analysis will be run from! You find more information about that in :ref:`Before Running <before_running>`.


    * -> `PACMAN/src/pacman/data/pacman_run_files/run_pacman.py <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/pacman_run_files/run_pacman.py>`_

    This python script runs PACMAN. If the user did not install PACMAN properly as explained in the :ref:`installation walkthrough <installation>` (not recommended) he or she might have to add a path to the sys.path.append pointing to the PACMAN directory:

    .. code-block:: python

  	  sys.path.append('/home/zieba/Desktop/Projects/Open_source/PACMAN/')

    The path in ``sys.path.append`` has to be changed to the location of PACMAN on the user's device.

    .. note:: The path should point to ``/PACMAN/`` and not to ``/PACMAN/src/pacman/``.


    * -> `PACMAN/src/pacman/data/pacman_run_files/obs_par.pcf <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/pacman_run_files/obs_par.pcf>`_

    The PACMAN control file (pcf): the user sets parameters related to the analysis. E.g., which plots should be saved, the path to the data directory and others.
    A thorough explanation of all the parameters in the pcf can be found on the :ref:`PCF page<pcf>`.


    * -> `PACMAN/src/pacman/data/pacman_run_files/fit_par.txt <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/pacman_run_files/fit_par.txt>`_

    The fit_par.txt file is used in Stage 30 to fit the light curve. The user defines here which fit parameters should be fixed, shared across visits, and sets other information like priors.




Nomenclature
''''''''''''''''''''''''''''''

An example for a directory structure:

::

    user
    └── Desktop
        └── Projects
            └── Observations
                └── Hubble
                    └── GJ1214_13021  (="run directory")
                        ├── pacman_run_files
                        │   ├── run_pacman.py
                        │   ├── obs_par.pcf
                        │   └── fit_par.txt
                        ├── stage00
                        ├── stage01
                        ├── stage02
                        ├── stage03
                        ├── stage10
                        ├── stage20
                        ├── stage21
                        └── stage30
            └── Data
                └── GJ1214_Hubble13021  (="data directory")
                    ├── ibxy06d0q_ima.fits
                    ├── ...
                    └── ibxy07ryq_ima.fits


* **run directory**:

  Example: ``/home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021``.

  This is the top-level directory for one PACMAN analysis. It contains the
  ``pacman_run_files`` directory and one output directory for each PACMAN stage.

* **pacman_run_files directory**:

  Example: ``/home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021/pacman_run_files``.

  This directory contains the files used to configure and start a PACMAN run:

  - ``run_pacman.py``
  - ``obs_par.pcf``
  - ``fit_par.txt``

  To run PACMAN, navigate into this directory and execute:

  .. code-block:: console

      python run_pacman.py

  The settings in ``obs_par.pcf`` and ``fit_par.txt`` are read from this directory.

* **stage directories**:

  PACMAN stores outputs in stage-specific directories:

  - ``stage00``
  - ``stage01``
  - ``stage02``
  - ``stage03``
  - ``stage10``
  - ``stage20``
  - ``stage21``
  - ``stage30``

  Each stage directory can contain one or more timestamped run directories. The
  name of each run directory is a combination of the stage name and the time at
  which the stage was started:

  ::

      s00_run_YYYY-MM-DD_HH-MM-SS

  For example:

  ::

      Stage00
      └── s00_run_2022-03-04_15-10-29

  When a later stage is run, PACMAN uses the most recent relevant run directory
  from the previous stage unless the user specifies a different input path.

* **data directory**:

  Example: ``/home/zieba/Desktop/Data/GJ1214_Hubble13021``.

  This directory contains the FITS files that will be reduced and analyzed.

* **pipeline directory**:

  Example: ``/home/zieba/Desktop/Projects/Open_Source/PACMAN/src``.

  This is the installed PACMAN source code directory. Users normally do not need to
  edit files in this directory during an analysis.

