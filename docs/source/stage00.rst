.. _stage00:

Stage 00
============

.. raw:: html

    <style> .blue {color:blue} </style>

.. topic:: TL;DR

    - Decide on a run directory
    - Copy pacman_script.py, fit_par.txt and obs_par.pcf to the run directory
    - Enter path to run directory and data directory into the pcf
    - Navigate to the run directory and execute the pacman_script.py file using the --s00 flag
    - Continue with s01


1) **Set up pcf**

    .. role:: blue
    :blue:`First, set up the location of your data directory (=datadir) and run directory (=rundir).`


    For me, these directories are located at the following places:

    - datadir = /home/zieba/Desktop/Data/GJ1214_Hubble13021
    - rundir = /home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021

    The rundir should contain the following three files:

    - **pacman_script.py**

    - **fit_par.txt**

    - **obs_par.pcf**

	These files can be either found in your package directory, on `GitHub <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files>`_
	or can be downloaded `here <https://downgit.github.io/#/home?url=https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files>`_.

	As mentioned in the Introduction, in this example going to analyze just two
	visits taken in the middle of the GO 13021 program (Dates (YYYY-MM-DD): 2013-03-13 and 2013-03-15).

    If your datadir contains a lot of visits, but you only want to analyze a subset you can use the ``which_visits`` parameter in the pcf.
    If the user downloaded all 15 visits in GO 13021, he or she can choose: ``which_visits = [5,6]`` in the pcf.


2) **Run PACMAN**

    You can give your analysis an "eventlabel" which might make it easier to identify it in the future.
    Here we just call our analysis: 'GJ1214_Hubble13021'.

	Now navigate to your rundir in your terminal and type:

	.. code-block:: console

		python pacman_script.py --s00 --eventlabel='GJ1214_Hubble13021'

	Here, --s00 means we are going to run Stage 00 and --eventlabel will be used in the naming of files and directories.
    The  --eventlabel argument is just needed when running --s00.

	When running s00, the first step creates a new subdirectory in rundir for the analysis which we will call the work directory (=workdir).

	It will be saved in the rundir and have a form like:
	``/home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021/run_2022-03-04_15-10-29_GJ1214_Hubble13021``

	The fit_par.txt and obs_par.pcf files which are in the run directory will be copied there.

	You should have a data structure similar to the following now:

::

	user
	└── Desktop
		└── Projects
			└── Observations
				└── Hubble
					└── GJ1214_13021
						├── pacman_script.py
						├── fit_par.txt
						├── obs_par.pcf
						└── run_2022-03-04_15-10-29_GJ1214_Hubble13021
							├── fit_par.txt
							└── obs_par.pcf



.. note::
    All next steps are going to use the pcf and fit_par which is located in the workdir (run_2022-03-04_15-10-29_GJ1214_Hubble13021) and not the pcf and fit_par in the rundir (GJ1214_13021)!!


3) Results

	After running Stage 00 you should get an output in the terminal similar to this one:

	.. code-block:: console

		    Starting s00
		    Found 1145 data file(s) ending in ima.fits
		    Reading in files and their headers: 100%|##########| 1145/1145 [00:03<00:00, 303.42it/s]
		    Determining orbit(s) and visit(s): 100%|##########| 1145/1145 [00:00<00:00, 261786.76it/s]
		    The user does not want to analyse every visit (which_visits != everything). The amount of files analyzed therefore reduced from 1145 to 150.
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
