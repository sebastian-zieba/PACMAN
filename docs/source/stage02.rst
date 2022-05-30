.. _stage02:

Stage 02
============

This stage performs a barycentric correction because the time stamps for the observations do not account for HST's orbital motion.

.. code-block:: console

	    Successfully reloaded meta file
	    Starting s02
	    Converting MJD to BJD: 100%|##########| 3/3 [00:02<00:00,  1.35it/s]
	    Writing t_bjd into filelist.txt
	    Saving Metadata
	    Finished s02


After the calculation has been performed, the user can check a newly generated plot also saved into "ancil/horizons".
Here we show the plot generated for the second of the two visits:

.. image:: media/s02/bjdcorr_horizons_results_v1.png

The axis are the distance of HST to the Solar System Barycenter in kilometers.
Horizons start and Horizons end show where our Horizon file starts and ends containing X,Y,Z information.
The black crosses in the plot show the times when HST actually observed. One can see that HST observed 4 orbits in this particular visit (which agrees with the ``filetable.txt`` from Stage 00).
One can also see the colored curve is a bit wiggley. This is in fact the rotation of HST around the earth.
The colored curve consists of a lot of points. Each one is an X,Y,Z position of HST downloaded from HORIZONS. The color coding denotes the time direction.

The ``filetable.txt`` is updated in this stage and contains a new column called ``t_bjd`` with the time of observations in BJD.
E.g. (only showing the first few lines):

.. include:: media/s02/filelist_updated.txt

   :literal:
