.. _stage00:

Stage 00
============

Let's start the analysis.
We are going to analyse just three visits taken in the middle of the program for simplicity.
Dates (YYYY-MM-DD): 2013-03-13,  2013-03-15,  2013-03-27
If the user has all 15 visits of the program in the data directory you can choose: ``which_visits = [5,6,7]`` in the pcf.

This step first creates a new directory for the analysis which will be used as the work directory ('workdir').
It will be saved in the rundir and have a form like:
``/home/zieba/Desktop/Projects/Open_source/PACMAN/run/run_2022-02-15_21-43-12_GJ1214_Hubble13021``

The fit_par.txt and obs_par.pcf files will be copied there.
After running Stage 00 you should get an output like this:


.. code-block:: console

	    Starting s00
	    Found 1145 data file(s) ending in ima.fits
	    Reading in files and their headers: 100%|##########| 1145/1145 [00:03<00:00, 303.42it/s]
	    Determining orbit(s) and visit(s): 100%|##########| 1145/1145 [00:00<00:00, 261786.76it/s]
	    The user does not want to analyse every visit (which_visits != everything). The amount of files analyzed therefore reduced from 1145 to 237.
	    Writing table into filelist.txt
	    Saving Metadata
	    Finished s00


You will also end up with a new file called ``filelist.txt``. It should look like this:

.. include:: media/s00/filelist.txt
   :literal:

It has the following columns:

* filenames

* instr: The specific filter or grism used in the obervation

* ivisit: Current visit of the observation

* iorbit: Current orbit of the observation

* t_mjd: Time in Modified Julian Date (MJD)

* t_visit: Time since the first exposure in the visit in minutes

* t_orbit: Time since the first exposure in the orbit in minutes

* scan: Scan direction:

  * 0: forward scan

  * 1: reverse scan

  * -1: not a spectrum but a direct image

* exp: exposure time in seconds
