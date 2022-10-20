.. _stage01:

Stage 01
============

.. topic:: TL;DR

    - Navigate to the run directory and execute the pacman_script.py file using the --s01 flag
    - Continue with s02

Next we download the locations of HST. This will be later used for the barycentric correction.

    .. warning:: This step needs an internet connection!

    .. note:: | At the beginning of every stage we read in again the pcf file located in the work directory.
              | This ensures that any user-made changes to the pcf file will be considered when running a new stage.
              | This means that the pcf file in the run directory is ONLY used in Stage 00. The same is true for the fit_par.txt file. So, after running Stage 00, PACMAN does not care anymore about the changes made to the pcf file and the fit_par file in the run directory.

Navigate to your rundir (where pacman_script.py is located) in your terminal and type:

.. code-block:: console

    python pacman_script.py --s01

The script (pacman_script.py) will assume that you want to continue the analysis started in the newly created work directory.
If you want to use a different work directory instead, you can use the --workdir='SOME_PATH' flag when running the script.


After running Stage 01 you should get an output like this:

.. code-block:: console

	    Successfully reloaded meta file
	    Starting s01
	    Retrieving Horizons file for every visit: 100%|##########| 2/2 [00:06<00:00,  2.04s/it]
	    Saving Metadata
	    Finished s01

We now accessed the `HORIZONS system <https://ssd.jpl.nasa.gov/horizons/>`_ by JPL and downloaded a file containing the positions of HST during the observations.
For that a new directory was created in the run directory called "ancil/horizons".
Two new .txt files where saved there; a Horizons file for each visit.
Each file contains the X, Y and Z position of HST relative to the solar system barycenter. The X,Y,Z positions of HST were downloaded for 5 minute intervals starting one hour before the first exposure in the observations and one hour after the observations.

For example, the first of the two horizon files should look like this (due to its length, we just display the first 100 lines):

.. include:: media/s01/horizons_results_v0_short.txt
   :literal:

The next Stage uses the information in these files to convert from MJD to BJD.

    .. note:: You might have noticed the output "Successfully reloaded meta file" at the beginning of the stage. This means that the pcf in the work directory is being read in again and any changes which have been made to the file between Stage 00 and Stage 01 will be considered. This reloading is being done before running every Stage but Stage 00.
