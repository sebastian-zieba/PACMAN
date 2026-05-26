.. _stage10:

Stage 10
============

.. topic:: Quick Summary

    - Navigate to ``pacman_run_files``, comment out Stage 03, uncomment Stage 10, and execute ``run_pacman.py``.
    - Continue with s20

1) **Stage summary**
This stage determines the position of the star in each direct image.

Let's look at an example using the two GJ1214 visits from earlier:

You can run Stage 10 first without giving any guess on where the star is located.
This will save a plot showing you the direct image which you can then use to refine your guess.

2) **Run PACMAN**

Navigate to ``pacman_run_files`` and open ``run_pacman.py``.
Comment out Stage 03 and uncomment Stage 10:

.. code-block:: python

    # meta = s03.run03(pcf_path=pcf_path)

    meta = s10.run10(pcf_path=pcf_path)



First run
''''''''''''''''''''''''

For demonstration purposes, in this first run I did not change the settings in the pcf file from a previous analysis of a different dataset:

| di_rmin    | 320
| di_rmax    | 360
| di_cmin    | 100
| di_cmax    | 150


Then run:

.. code-block:: console

    python run_pacman.py


The terminal should give you something like this:

.. code-block:: console

	Starting s10
	Using Stage 03 input directory: /Users/sebastianzieba/Desktop/Projects/Observations/Hubble/GJ1214_13021_2026/stage03/s03_run_2026-05-26_12-51-37
	Location of the new Stage 10 run directory: /Users/sebastianzieba/Desktop/Projects/Observations/Hubble/GJ1214_13021_2026/stage10/s10_run_2026-05-26_12-57-19
	Determining Source Positions for Direct Images:   0%|        | 0/8 [00:00<?, ?it/s]
	Your guess for di_rmax is outside of the image.
	Determining Source Positions for Direct Images:  12%|#       | 1/8 [00:00<00:02,  2.34it/s]
	Your guess for di_rmax is outside of the image.
	Determining Source Positions for Direct Images:  25%|##      | 2/8 [00:00<00:01,  3.03it/s]
	Your guess for di_rmax is outside of the image.
	Determining Source Positions for Direct Images:  38%|###     | 3/8 [00:00<00:01,  3.41it/s]
	Your guess for di_rmax is outside of the image.
	Determining Source Positions for Direct Images:  50%|####    | 4/8 [00:01<00:01,  3.62it/s]
	Your guess for di_rmax is outside of the image.
	Determining Source Positions for Direct Images:  62%|#####   | 5/8 [00:01<00:00,  3.74it/s]
	Your guess for di_rmax is outside of the image.
	Determining Source Positions for Direct Images:  75%|######  | 6/8 [00:01<00:00,  3.83it/s]
	Your guess for di_rmax is outside of the image.
	Determining Source Positions for Direct Images:  88%|####### | 7/8 [00:01<00:00,  3.89it/s]
	Your guess for di_rmax is outside of the image.
	Determining Source Positions for Direct Images: 100%|########| 8/8 [00:02<00:00,  3.64it/s]
	There is one DI per orbit.
	Saving Metadata
	Finished s10 

We see that we got error messages that our guesses were outside of the image.
To determine where the star actually is, this run saved the plots into the current Stage 10 workdir under:

``stage10/s10_run_*/figs/images``

You will end up with plots like this one:

.. image:: media/s10/quick_di0_wrong.png

You can see that our "initial guess" (the red box in the plot) was off by a lot from the star. Let's try again with a better guess, now that we know where the star is.


Second run
''''''''''''''''''''''''

By estimating by eye we can tell the star is approximately at row = 140 and col = 30. So our new guess is:

| di_rmin   | 120
| di_rmax   | 160
| di_cmin   | 5
| di_cmax   | 50


The terminal should give you something like this:

.. code-block:: console

	Starting s10
	Using Stage 03 input directory: /Users/sebastianzieba/Desktop/Projects/Observations/Hubble/GJ1214_13021_2026/stage03/s03_run_2026-05-26_12-51-37
	Location of the new Stage 10 run directory: /Users/sebastianzieba/Desktop/Projects/Observations/Hubble/GJ1214_13021_2026/stage10/s10_run_2026-05-26_13-00-32
	Determining Source Positions for Direct Images: 100%|########| 8/8 [00:05<00:00,  1.38it/s]
	There is one DI per orbit.
	Saving Metadata
	Finished s10 


You will end up with plots like this one:

.. image:: media/s10/quick_di0.png

This time we have the star in our cutout box!
A second plot shows you the best fit to the star by using a 2D gaussian.

.. image:: media/s10/di_0.png

The positions of the star in physical pixels are saved into:

``stage10/s10_run_*/xrefyref.txt``

.. include:: media/s10/xrefyref.txt
   :literal:

| 1. column is the time when the direct image was taken in BJD
| 2. column is the visit number
| 3. column is the orbit number
| 4. column is the cumulative orbit number
| 5. column is the row position of the star in physical pixels
| 6. column is the column position of the star in physical pixels

You might notice that the fit resulted in a target location of approximately row = 513 and col = 403 but in the plots the target is closer to row = 140 and col = 30 (estimating by eye).
This is because the table accounted for the offset in X and Y to subarray start. This information can be found in the header:

| LTV1    =  -374.0 / offset in X to subsection start
| LTV2    =  -374.0 / offset in Y to subsection start

Also note that you have two s10 runs in your Stage 10 directory. You can delete the first one if you want to save space. PACMAN will not use the first run for the next stages, it will automatically use the most recent one. So as long as you have the second run in there, you are good to go for the next stages.
