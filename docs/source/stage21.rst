.. _stage21:

Stage 21
============

.. topic:: Summary

    - Navigate to the run directory and execute the run_pacman.py file using the --s21 flag
    - Continue with s30 to fit the white and spectroscopic light curves

This stage bins the data into spectroscopic light curves. 

PACMAN automatically uses the most recent Stage 20 run.

The user can set the wavelength bins in the pcf by giving a wavelength range and a number of bins, or alternatively provide user-defined bin edges.


When running S21 the user should get an output similar to this one:

.. code-block:: console

    Starting s21
    Using Stage 20 input directory: ...
    Location of the new Stage 21 run directory: ...

    Number of bins: 11
    chosen bin edges: [11370.         11830.90909091 12291.81818182 12752.72727273
     13213.63636364 13674.54545455 14135.45454545 14596.36363636
     15057.27272727 15518.18181818 15979.09090909 16440.        ]
    ***************** Looping over Bins: 100%|########################################| 11/11 [00:01<00:00,  9.15it/s]
    Saved light curve(s) in stage21/s21_run_*/extracted_sp
    Saving Wavelength bin file
    Saving Metadata
    Finished s21


Below is a plot of a 1D spectrum with the bin edges from a user-chosen binning.

.. image:: media/s21/spec_bins11.png

The spectroscopic light curves are saved in:

``stage21/s21_run_*/extracted_sp``

For each wavelength bin, PACMAN saves:

- a spectroscopic light curve (``speclc*.txt``)
- light-curve plots with uncertainties
- individual visit light curves (if multiple visits exist)

The ``extracted_sp`` directory also contains:

- ``wvl_table.dat``:
  wavelength-bin information including:
  
  - bin number
  - central wavelength
  - half-width
  - lower wavelength edge
  - upper wavelength edge

- ``figs/spec_binsXX.png``:
  overview plot showing the wavelength bins on the spectrum

The wavelength-bin figures are saved in:

``stage21/s21_run_*/extracted_sp/figs``

The stucture will be something like this:

::

  stage21/s21_run_YYYY-MM-DD_HH-MM-SS/
  ├── extracted_sp/
  │   ├── speclc1.158.txt
  │   ├── speclc1.204.txt
  │   ├── ...
  │   └── wvl_table.dat
  └── figs/
      └── ...
      