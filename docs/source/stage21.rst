.. _stage21:

Stage 21
============

This stage bins the data into spectroscopic light curves. 

PACMAN can either use the most recent S20 run for the S21 step or specify a path to a previous run.

The user can set the wavelength bins in the pcf by giving a wavelength range and a number of bins, or alternatively provide user-defined bin edges.


When running S21 the user should get an output similar to this one:

.. code-block:: console

	    Successfully reloaded meta file
	    Starting s21

	    Number of bins: 12
	    chosen bin edges: [11250.  11687.5 12125.  12562.5 13000.  13437.5 13875.  14312.5 14750.
	    15187.5 15625.  16062.5 16500. ]
	    Chosen directory with the spectroscopic flux file: 2022-02-15_22-08-24
	    ***************** Looping over Bins: 100%|##########| 12/12 [00:05<00:00,  2.25it/s]
	    saved light curve(s) in run_2022-02-15_21-52-46_GJ1214_Hubble13021//extracted_sp/bins12_2022-02-15_22-14-31
	    Finished s21

Below is a plot of a 1D spectrum with the bin edges from a user-chosen binning.

.. image:: media/s21/spec_bins12.png
