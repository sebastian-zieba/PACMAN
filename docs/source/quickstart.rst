.. _quickstart:

Quickstart
============

1. Installation and requirements
-------------------------------------------

See installation.rst


2. Download data
-------------------------------------------

ima files!


3. Setup the PACMAN control file (.pcf)
-------------------------------------------

See pcf.rst

4. Run PACMAN
-------------------------------------------

PACMAN is separated into different stages. Here a quick summary before we look into each of them ińto more detail by
applying PACMAN on real observations:

- Stage 00: Every observation in the data directory is being read in an import information from the header are saved into a new file.

- Stage 01: To prepare for the barycentric correction, we download the location of HST during the observations

- Stage 02: We perform the barycentric correction

- Stage 03: We load in a stellar spectrum. The product of the stellar spectrum and the bandpass of the Grism (G102 or G141) will be used as a reference spectrum in the wavelength calibration.

- Stage 10: The direct images are being read in and the position of the star is being determined.

- Stage 20: We extract the spectra

- Stage 21: TODO! Generate limb darkening parameters

- Stage 22: Bin the spectra to create light curves

- Stage 30: Fit a model to the light curves


Introduction
:::::::::::::::::::::::::::::::::::::::::

Let's apply PACMAN on real observations of GJ1214. The planet was observed in HST GO13021 in 15 visits.
Let's just look at the last two visits (so which_visits = [13,14]) for simplicity.

Nomenclature:
:::::::::::::::::::::::::::::::::::::::::
* run directory:
/home/zieba/Desktop/Projects/Open_source/PACMAN/run

* work directory:
/home/zieba/Desktop/Projects/Open_source/PACMAN/run/run_2022-01-19_16-46-19_GJ1214_Hubble13021

* data directory:
/home/zieba/Desktop/Data/GJ1214_Hubble13021

Stage 00
:::::::::::::::::::::::::::::::::::::::::

This step first creates a new directory for the analysis which will be used as the work directory ('workdir').
It will be saved in the rundir and have a form like:
/home/zieba/Desktop/Projects/Open_source/PACMAN/run/run_2022-01-19_16-46-19_GJ1214_Hubble13021

The fit_par.txt and obs_par.pcf files will be copied there.
After running Stage 00 you should get an output like this:

| Starting s00
| Found 1145 data file(s) ending in ima.fits
| Reading in files and their header information: 100%|██████████| 1145/1145 [00:03<00:00, 368.21it/s]
| Determining Orbit and Visit: 100%|██████████| 1145/1145 [00:00<00:00, 311606.42it/s]
| The user does not want to analyse every visit (which_visits != everything). The amount of files analyzed therefore reduced from 1145 to 158.
| Writing table into filelist.txt
| Saving Metadata
| Finished s00

You will also end up with a new file called filelist.txt. It should look like this:

.. include:: media/filelist.txt
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


Stage 01
:::::::::::::::::::::::::::::::::::::::::


Next we download the locations of HST. This will be later used for the barycentric correction.

    .. warning:: This step needs an internet connection!

    .. note:: | At the beginning of every stage we read in again the pcf file located in the work directory (e.g. /home/zieba/Desktop/Projects/Open_source/PACMAN/run/run_2022-01-19_16-46-19_GJ1214_Hubble13021).
              | This ensures that any user-made changes to the pcf file will be considered when running a new stage.
              | This means that the pcf file in the run directory (e.g., /home/zieba/Desktop/Projects/Open_source/PACMAN/run) is ONLY used in Stage 00. The same is true for the fit_par.txt file. So, after running Stage 00, PACMAN does not care anymore about the changes made to the pcf file and the fit_par file in the run directory.

After running Stage 01 you should get an output like this:

.. code-block:: console

	    Successfully reloaded meta file
	    Starting s01
	    Retrieving Horizons file for every visit: 100%|██████████| 2/2 [00:02<00:00,  1.03s/it]
	    Saving Metadata
	    Finished s01

We now accessed the `HORIZONS system <https://ssd.jpl.nasa.gov/horizons/>`_ by JPL and downloaded a file containing the positions of HST during the observations.
For that a new directory was created in the run directory called "ancil/horizons".
Two new .txt files where saved there; a Horizons file for each visit.
Each file contains the X, Y and Z position of HST relative to the solar system barycenter. The X,Y,Z positions of HST were downloaded for 5 minute intervals starting one hour before the first exposure in the observations and one hour after the observations.

For example, the second file should look like this (due to its length, we just display the first 100 lines of it):

.. include:: media/horizons_results_v1_short.txt
   :literal:

The next Stage uses the information in these files to convert from MJD to BJD.

Stage 02
:::::::::::::::::::::::::::::::::::::::::

This Stage has to perform a barycentric correction because the header only contains MJD.

.. code-block:: console

	    Successfully reloaded meta file
	    Starting s02
	    Converting MJD to BJD: 100%|██████████| 2/2 [00:01<00:00,  1.27it/s]
	    Writing t_bjd into filelist.txt
	    Saving Metadata
	    Finished s02

After the calculation has been performed, the user can check a newly generated plot also saved into "ancil/horizons".
Here we show the plot generated for the second of the two visits:

.. image:: media/bjdcorr_horizons_results_v1.png

The axis are the distance of HST to the Solar System Barycenter in kilometers.
Horizons start and Horizons end show where our Horizon file starts and ends containing X,Y,Z information.
The black crosses in the plot show the times when HST actually observed. One can see that HST observed 4 orbits in this particular visit (which agrees with the filetable.txt from Stage 00).
One can also see the colored curve is a bit wiggley. This is in fact the rotation of HST around the earth.
The colored curve consists out of a lot of points. Every one of them are a X,Y,Z position of HST downloaded from HORIZONS. The color coding denotes the time direction.

Stage 03
:::::::::::::::::::::::::::::::::::::::::

This Stage starts by downloading a stellar model.

    .. note:: | PACMAN currently offers the following stellar models:
              | - THE 1993 KURUCZ STELLAR ATMOSPHERES ATLAS (k93models)
              | - THE STELLAR ATMOSPHERE MODELS BY CASTELLI AND KURUCZ 2004 (ck04models)
              | - THE PHOENIX MODELS BY FRANCE ALLARD AND COLLABORATORS (phoenix)
              | - blackbody spectrum
              | The stellar models (exluding the blackbody one) are retrieved from https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/

We use the stellar parameters published in `Cloutier et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021AJ....162..174C/abstract>`_.

| Teff   =  3250
| logg   =  5.026
| MH     =  0.29

.. code-block:: console

	    Successfully reloaded meta file
	    Starting s03
	    Using g141 grism.
	    Using phoenix model.

	    Possible metallicities: [ 0.  -0.5 -1.  -1.5 -2.  -2.5 -3.  -3.5 -4.   0.5  0.3]
	    Using input metallicity of -0.5.

	    Possible effective temperatures: [10000. 10200. 10400. 10600. 10800. 11000. 11200. 11400. 11600. 11800.
	     12000. 12500. 13000. 13500. 14000. 14500. 15000. 15500. 16000. 16500.
	     17000. 17500. 18000. 18500. 19000. 19500.  2000. 20000.  2100. 21000.
	      2200. 22000.  2300. 23000.  2400. 24000.  2500. 25000.  2600. 26000.
	      2700. 27000.  2800. 28000.  2900. 29000.  3000. 30000.  3100. 31000.
	      3200. 32000.  3300. 33000.  3400. 34000.  3500. 35000.  3600. 36000.
	      3700. 37000.  3800. 38000.  3900. 39000.  4000. 40000.  4100. 41000.
	      4200. 42000.  4300. 43000.  4400. 44000.  4500. 45000.  4600. 46000.
	      4700. 47000.  4800. 48000.  4900. 49000.  5000. 50000.  5100. 51000.
	      5200. 52000.  5300. 53000.  5400. 54000.  5500. 55000.  5600. 56000.
	      5700. 57000.  5800. 58000.  5900. 59000.  6000. 60000.  6100. 61000.
	      6200. 62000.  6300. 63000.  6400. 64000.  6500. 65000.  6600. 66000.
	      6700. 67000.  6800. 68000.  6900. 69000.  7000. 70000.  7200.  7400.
	      7600.  7800.  8000.  8200.  8400.  8600.  8800.  9000.  9200.  9400.
	      9600.  9800.]
	    For input effective temperature 3412, closest temperature is 3400.0.

	    Was the stellar model fits file already downloaded?: False

		          + Downloading file phoenixm05_3400.fits from https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/phoenix/phoenixm05/phoenixm05_3400.fits.
	    --2022-01-20 02:23:51--  https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/phoenix/phoenixm05/phoenixm05_3400.fits
	    Resolving archive.stsci.edu (archive.stsci.edu)... 130.167.201.60
	    Connecting to archive.stsci.edu (archive.stsci.edu)|130.167.201.60|:443... connected.
	    HTTP request sent, awaiting response... 200 OK
	    Length: 532800 (520K) [image/fits]
	    Saving to: ‘phoenixm05_3400.fits’

	         0K .......... .......... .......... .......... ..........  9%  170K 3s
	        50K .......... .......... .......... .......... .......... 19%  374K 2s
	       100K .......... .......... .......... .......... .......... 28% 1.27M 1s
	       150K .......... .......... .......... .......... .......... 38% 1.51M 1s
	       200K .......... .......... .......... .......... .......... 48%  376K 1s
	       250K .......... .......... .......... .......... .......... 57% 5.90M 0s
	       300K .......... .......... .......... .......... .......... 67% 1.86M 0s
	       350K .......... .......... .......... .......... .......... 76%  344K 0s
	       400K .......... .......... .......... .......... .......... 86% 1.25M 0s
	       450K .......... .......... .......... .......... .......... 96% 1003K 0s
	       500K .......... ..........                                 100% 16.8M=0.9s

	    2022-01-20 02:23:53 (577 KB/s) - ‘phoenixm05_3400.fits’ saved [532800/532800]

	    Possible logg: [0.  0.5 1.  1.5 2.  2.5 3.  3.5 4.  4.5 5.  5.5]
	    For input logg 4.94, closest logg is 5.0.

	    Saving Metadata
	    Finished s03


.. image:: media/refspec.png

