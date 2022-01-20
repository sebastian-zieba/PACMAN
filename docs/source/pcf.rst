.. _pcf:

PACMAN Control File (.pcf)
============================

To run the different Stages of ``PACMAN``, the pipeline requires control files (.pcf) where Stage-specific parameters are defined (e.g. aperture size, path of the data, etc.).

In the following, we look at the contents of the pcf.



.. include:: media/obs_par.pcf
   :literal:


Stage 00
---------------------------------------------------------

rundir
''''''''''''''''''''''''''''''''''''''''''''
Example: ``/home/zieba/Desktop/Projects/Open_source/PACMAN/run``
The directory where you want PACMAN to run and save data to.
If you downloaded or cloned the GitHub repository it includes a run directory you can use for the runs.
It has to include three files:
 - run_pacman.py: The run script
 - obs_par.pcf: The pcf file
 - fit_par.txt: The fit_par file with the fit parameters (if you don't want to fit the data (= Stage 30) you dont have to care about the contents of this file)

.. warning::

	Your path is not allowed to have any spaces in it. E.g., ``/home/USER/run 1`` is not a valid path.

datadir
''''''''''''''''''''''''''''''''''''''''''''
Example: ``/home/zieba/Desktop/Data/GJ1214_Hubble13021``
The location of your data

.. warning::

	Your path is not allowed to have any spaces in it. E.g., ``/home/USER/data GJ1214`` is not a valid path.

suffix
''''''''''''''''''''''''''''''''''''''''''''
Only possible extended supported currently: ``ima``
From the `WFC3 data handbook (Types of WFC3 Files) <https://hst-docs.stsci.edu/wfc3dhb/chapter-2-wfc3-data-structure/2-1-types-of-wfc3-files>`_: "For the IR detector, an intermediate MultiAccum (ima) file is the result after all calibrations are applied (dark subtraction, linearity correction, flat fielding, etc.) to all of the individual readouts of the IR exposure."

which_visits
''''''''''''''''''''''''''''''''''''''''''''
| Example: ``[0,2]``
| Example: ``everything``

If your ``datadir`` contains several visits, you can select which ones you want to analyse.
If you are interested in all visits in ``datadir`` use ``everything`` here.

E.g., HST GO 13021 contains 15 visits. If you have all 15 visits in your ``datadir`` but only want to analyse the last two visits for now, you would enter ``[13,14]`` here.
If your ``datadir`` only contained these two visits (and not the previous 13 visits before it), you can either write ``everything`` or ``[0,1]``.



Stage 02
---------------------------------------------------------

save_barycorr_plot/show_barycorr_plot
''''''''''''''''''''''''''''''''''''''''''''

Saves or shows a plot with the downloaded X,Y,Z positions of HST from the `HORIZONS system <https://ssd.jpl.nasa.gov/horizons/>`_ by JPL during the observations.

.. image:: media/bjdcorr_horizons_results_v1.png




Stage 03
---------------------------------------------------------

Teff, logg, MH
''''''''''''''''''''''''''''''''''''''''''''
| Example: ``Teff    3250``
| Example: ``logg    5.026``
| Example: ``MH      0.29``

effective Temperature (Teff), surface gravity (logg) and metallicity (MH) of the star.

sm
''''''''''''''''''''''''''''''''''''''''''''
Example: ``phoenix``

The stellar model one wants to use. It will be then multiplied with the bandpass (G102 or G141) to create a reference spectrum for the wavelength calibration of the spectra.

Options:

| PACMAN currently offers the following stellar models:
| - ``k93models``: THE 1993 KURUCZ STELLAR ATMOSPHERES ATLAS
| - ``ck04models``: THE STELLAR ATMOSPHERE MODELS BY CASTELLI AND KURUCZ 2004
| - ``phoenix``: THE PHOENIX MODELS BY FRANCE ALLARD AND COLLABORATORS
| - ``blackbody``: A blackbody spectrum usings Planck's law

The stellar models (exluding the blackbody one) are retrieved from https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/


save_refspec_plot/show_refspec_plot
''''''''''''''''''''''''''''''''''''''''''''

.. image:: media/refspec.png



Stage 10
---------------------------------------------------------

di_rmin, di_rmax, di_cmin, di_cmax
''''''''''''''''''''''''''''''''''''''''''''


save_image_plot/show_image_plot
''''''''''''''''''''''''''''''''''''''''''''


