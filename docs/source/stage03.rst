.. _stage03:

Stage 03
============

.. topic:: TL;DR

    - Navigate to the run directory and execute the pacman_script.py file using the --s03 flag
    - Continue with s10

This Stage starts by downloading a stellar model or calculating one in case of the black body spectrum.

    .. note:: | PACMAN currently offers the following stellar models:
              | - THE 1993 KURUCZ STELLAR ATMOSPHERES ATLAS (``k93models``)
              | - THE STELLAR ATMOSPHERE MODELS BY CASTELLI AND KURUCZ 2004 (``ck04models``)
              | - THE PHOENIX MODELS BY FRANCE ALLARD AND COLLABORATORS (``phoenix``)
              | - ``blackbody`` spectrum
              | The stellar models (exluding the blackbody) are retrieved from https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/

    .. warning:: If the user decides on any stellar model which has to be downloaded (``k93models``, ``ck04models`` or ``phoenix``), then internet connection is required! The blackbody model does not require a connection however.

The user sets Teff, logg and MH in the pcf file. For ``k93models``, ``ck04models`` and ``phoenix`` the closest available metallicity,
then the closest available effective temperature and finally the closest available surface gravity is searched and then used.
To see the workflow see the code block below.

Let's look at an example for GJ 1214:
We use the stellar parameters published in `Cloutier et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021AJ....162..174C/abstract>`_.

| Teff   =  3250
| logg   =  5.026
| MH     =  0.29

Let's use an ``k93models`` model.

.. code-block:: console

	    Successfully reloaded meta file
	    Starting s03
	    Using k93models model.

	    Possible metallicities: [ 1.   0.5  0.3  0.2  0.1  0.  -0.1 -0.2 -0.3 -0.5 -1.  -1.5 -2.  -2.5
	     -3.  -3.5 -4.  -4.5 -5. ]
	    For input metallicity 0.29, closest metallicity is 0.3.

	    Possible effective temperatures: [10000. 10500. 11000. 11500. 12000. 12500. 13000. 14000. 15000. 16000.
	     17000. 18000. 19000. 20000. 21000. 22000. 23000. 24000. 25000. 26000.
	     27000. 28000. 29000. 30000. 31000. 32000. 33000. 34000.  3500. 35000.
	      3750. 37500.  4000. 40000.  4250. 42500.  4500. 45000.  4750. 47500.
	      5000. 50000.  5250.  5500.  5750.  6000.  6250.  6500.  6750.  7000.
	      7250.  7500.  7750.  8000.  8250.  8500.  8750.  9000.  9250.  9500.
	      9750.]
	    For input effective temperature 3250, closest temperature is 3500.0.

	    Was the stellar model fits file already downloaded?: False

		          + Downloading file kp03_3500.fits from https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/k93models/kp03/kp03_3500.fits.
	    --2022-01-25 19:13:30--  https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/k93models/kp03/kp03_3500.fits
	    Resolving archive.stsci.edu (archive.stsci.edu)... 130.167.201.60
	    Connecting to archive.stsci.edu (archive.stsci.edu)|130.167.201.60|:443... connected.
	    HTTP request sent, awaiting response... 200 OK
	    Length: 69120 (68K) [image/fits]
	    Saving to: ‘kp03_3500.fits’

	    kp03_3500.fits         100%[=========================>]  67.50K   323KB/s    in 0.2s

	    2022-01-25 19:13:30 (323 KB/s) - ‘kp03_3500.fits’ saved [69120/69120]

	    Possible logg: [0.  0.5 1.  1.5 2.  2.5 3.  3.5 4.  4.5 5. ]
	    For input logg 5.026, closest logg is 5.0.

	    Using g141 grism.
	    Saving Metadata
	    Finished s03



.. note::

	The stellar models will be saved into the run directory. E.g., ``/home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021/run_2022-03-04_15-10-29_GJ1214_Hubble13021/ancil/stellar_models/k93models/kp03_3500.fits``.
	If the file already exists then it will not be downloaded again.


After downloading or calculating the stellar spectrum, it is multiplied by the grism throughput (if the observations used G102 or G141 is recognized automatically).
The bandpass files are stored in the pipeline directory (e.g., ``/home/zieba/Desktop/Projects/Open_source/PACMAN/pacman/ancil/bandpass``).
The final plot of this stage shows the stellar spectrum, the bandpass and the product of these two. This product is used as the reference spectrum for wavelength calibration.

.. image:: media/s03/refspec.png
