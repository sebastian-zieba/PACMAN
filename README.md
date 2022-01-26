![](https://github.com/sebastian-zieba/PACMAN/blob/master/docs/source/media/Pacman_V1.gif)

# PACMAN

**Welcome to PACMAN**

ALERT: Project PACMAN is currently under heavy development. Use at your own risk.

## Installation

TODO: Add short installation guide

### With Git/GitHub

1. You can install PACMAN directly from source on [GitHub](https://github.com/sebastian-zieba/PACMAN) in one of two ways:
	- On the [GitHub website](https://github.com/sebastian-zieba/PACMAN), click on **Code** and **Download ZIP** followed by unpacking the distribution by opening up a terminal and typing:

		```bash
		unzip PACMAN-master.zip
		```

	- OR, clone the repository using ``git`` by typing:

		```bash
		git clone https://github.com/sebastian-zieba/PACMAN.git
		```

2. Navigate into the newly created directory and **install** PACMAN by running ``setup.py``:

	```bash
	python setup.py install
	```

3. Install additional **requirements** for the package by typing:

	```bash
	pip install -r requirements.txt
	```


## Documentation

Check out the docs [here](https://pacmandocs.readthedocs.io/en/latest/)




## PACMAN Steps:



- ./wfc3_reduction/reduction/s00_table.py
  - Creates a new run directory with the following form, eg.: ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/
  - Copy and pastes the control file (obs_par.ecf) and the fit parameters file (fit_par_new2.txt) into the new directory
  - Reads in all fits files and creates a table which will be saved in filelist.txt. The following information will be listed in this table: filenames filter/grism ivisit iorbit t_mjd t_visit t_orbit scan exp

Meta data will be saved after every step in a file call something like ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/WFC3_L-98-59_Hubble15856_Meta_Save.dat
Meta data includes information like fits headers containing observation specific parameters.  
#
- ./wfc3_reduction/reduction/s01_horizons.py
  - Retrieves vector data of Hubble from JPL's HORIZONS system on https://ssd.jpl.nasa.gov/horizons_batch.cgi (see Web interface on https://ssd.jpl.nasa.gov/horizons.cgi)
  - txt file with HST positions in space will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/ancil/horizons
#
- ./wfc3_reduction/reduction/s02_barycorr.py
  - performs barycorr based on the t_mjd in filelist.txt. Adds another column called t_bjd
  - Plots will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/ancil/horizons
#
- ./wfc3_reduction/reduction/s03_refspectra.py
  - Downloads the bandpass of G141 or G102 and a stellar spectrum using the python package pysynphot. The product of these two will be used as a reference spectrum and for the wavelength calibration.
  - Plots and data will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/ancil/bandpass
#
- ./wfc3_reduction/reduction/s10_direct_images.py
  - The direct images will be read in and the centroid position of the star will be determined.
  - This information will be saved in a new file called ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/xrefyref.txt
  - A plot of the direct image will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/figs/images
#
- ./wfc3_reduction/reduction/s20_extract.py
  - Creates a new directory called ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/extracted_lc/2021-08-03_12:22/
  - The white lightcurve will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/extracted_lc/2021-08-03_12:22/lc_white.txt
  - The spectroscopic lightcurve data will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/extracted_lc/2021-08-03_12:22/lc_spec.txt (the actual binning happens in step s22)
  - Several plots will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/figs/images
#
- ./wfc3_reduction/reduction/s21_ld_inputmaker.py OPTIONAL
  - Doesnt really completely work yet! Will create limb darkening parameters if the user doesnt want to fit for them but fix them to theoretical values.
#
- ./wfc3_reduction/reduction/s22_bin_spectroscopic_lc.py
  - Creates a new directory called ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/extracted_sp/bins12_2021-08-09_04-36-06/ and saves the binned spectroscopic lcs there. 
#
- ./wfc3_reduction/reduction/s30_run.py
  - Fits the spectroscopic or white lcs. Set fit_par_new2.txt for that!

## PACMAN control file:

This file located in ./run has to be set before running any step! Let's have a look at the different parameters here:

Some short explanations:

### 00
- path     /home/zieba/Desktop/Data/L-98-59_Hubble15856/ #location of fits files
- suffix             ima 
- which_visits  [0] # just use the 0th visit

### 10
- di_rmin  270 #estimated position of the star in the direct image. that cutout will be used to determine the centroid position then
- di_rmax  290
- di_cmin  245
- di_cmax  265

### 20
- window   10: plus and minus the size of aperture in pixels of the row with the highest gradient in flux
- diagnostics        False
- background_thld  1000 #background threshold in counts. will just take the median below this value to determine the background flux
- opt_extract  True # want to use optimal extraction? If not the flux in the box aperture will be simply added up
- sig_cut  15  #optimal extraction, for cosmic rays etc
- nsmooth  9         #optimal extraction, created smoothed spatial profile, medial smoothing filter
- rdnoise  22.0 #constant for HST
- output   True
- correct_wave_shift   True

### 22
- wvl_min  1.125 #start of wavenlegth range to consider
- wvl_max  1.65 #end of wavenlegth range to consider
- wvl_bins   [6] #6 bins. user can also enter an array here like [6,12,18]. This will create three different directories with the different spec lcs then 

##30
- toffset  2458366 #time offset which will be subtracted from the pobservations to prevent a floating point precision problem

- run_verbose         True
- run_output          False
- run_show_plot       True
- run_mcmc            True  #Do MCMC with emcee
- run_nested          False #Do NS with dynesty
- run_lsq             True  #Do least square (with MPFit) should always be true
- run_plot_raw_data   True  #Create a plot with the raw light curve
- run_fit_white       True  #SET THIS TRUE IF YOU WANT TO DO A WHITE LC FIT. IF YOU WANT TO FIT THE SPECTR LCS SET FALSE 
- run_divide_white    False # i usually set this false. uses the divide white method. not tested

- run_files           ['/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-16_22-31-42_L-98-59_Hubble15856/extracted_sp/bins12_2021-08-17_23-34-29'] # example for spec fit
- run_files           ['/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-21_02-00-29_L-98-59_Hubble15856/extracted_lc/2021-08-21_02:01/lc_white.txt'] # example for white fit
- run_myfuncs         ['constant','upstream_downstream','model_ramp','polynomial1','transit'] #['constant', 'upstream_downstream', 'ackbar', 'polynomial1', 'transit']# ['constant', 'upstream_downstream', 'model_ramp', 'polynomial2', 'transit']


- save_allan_plot    False

####mcmc
- run_nsteps          12000
- run_nwalkers	     75
- run_nburn           6000

####Ns
- run_dlogz     25000
- run_nlive     100


- run_clipsigma       4 #will remove outliers iteratively and fit the lc. in this case 4 sigma outliers. it will iterate a max of 4 times. i believe that works just for the white lc. might be good to set both to 0 for now
- run_clipiters       4


