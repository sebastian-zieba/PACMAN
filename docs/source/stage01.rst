.. _stage01:

Stage 01
============

.. topic:: Quick Summary

    - Navigate to ``pacman_run_files``, only have Stage 01 uncommented and execute run_pacman.py.
    - Continue with s02

1) **Stage summary**
Next we download the locations of HST. This will be later used for the barycentric correction.

    .. warning:: This step needs an internet connection!

2) **Run PACMAN**
Navigate to ``pacman_run_files`` and open ``run_pacman.py``. Comment out Stage 00 and uncomment Stage 01:

.. code-block:: python

    # meta = s00.run00(pcf_path=pcf_path)

    meta = s01.run01(pcf_path=pcf_path)

Then run:

.. code-block:: console

    python run_pacman.py

3) **What happens?**
PACMAN reads the most recent ``stage00/s00_run_*`` directory and creates a new timestamped Stage 01 workdir:

``stage01/s01_run_YYYY-MM-DD_HH-MM-SS``

The HORIZONS files are saved in:

``stage01/s01_run_*/ancil/horizons``

After running Stage 01 you should get an output like this:

.. code-block:: console

    Starting s01
    Using Stage 00 input directory: /Users/sebastianzieba/Desktop/Projects/Observations/Hubble/GJ1214_13021_2026/stage00/s00_run_2026-05-26_12-33-34
    Location of the new Stage 01 run directory: /Users/sebastianzieba/Desktop/Projects/Observations/Hubble/GJ1214_13021_2026/stage01/s01_run_2026-05-26_12-41-21
    Retrieving Horizons file for every visit: 100%|##########| 2/2 [00:01<00:00,  1.95it/s]
    Saving Metadata
    Finished s01 

We now accessed the `HORIZONS system <https://ssd.jpl.nasa.gov/horizons/>`_ by JPL and downloaded a file containing the positions of HST during the observations.
The files are saved in the current Stage 01 workdir under ``ancil/horizons``.
Two new .txt files were saved there; a Horizons file for each visit.
Each file contains the X, Y and Z position of HST relative to the solar system barycenter. The X,Y,Z positions of HST were downloaded for 5 minute intervals starting one hour before the first exposure in the observations and one hour after the observations.

For example, the first of the two horizon files should look like this (due to its length, we just display the first 100 lines):

.. include:: media/s01/horizons_results_v0_short.txt
   :literal:

You should now have a directory structure similar to this:

TEST

::

    user
    └── Desktop
        └── Projects
            └── Observations
                └── Hubble
                    └── GJ1214_13021
                        ├── pacman_run_files
                        │   ├── run_pacman.py
                        │   ├── fit_par.txt
                        │   └── obs_par.pcf
                        ├── stage00
                        │   └── s00_run_2022-03-04_15-10-29
                        │       ├── fit_par.txt
                        │       ├── obs_par.pcf
                        │       ├── s00.log
                        │       ├── WFC3_Meta_Save.dat
                        │       └── filelist.txt
                        └── stage01
                            └── s01_run_2022-03-04_15-10-29
                                ├── fit_par.txt
                                ├── obs_par.pcf
                                ├── s01.log
                                ├── WFC3_Meta_Save.dat
                                └── filelist.txt
                                └── ancil
                                    └── horizons
                                        ├── horizons_results_v0.txt
                                        └── horizons_results_v1.txt

The next Stage uses the information in these files to convert from MJD to BJD.

.. include:: media/s01/horizons_results_v0_short.txt
   :literal:
