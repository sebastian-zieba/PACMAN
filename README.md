![](https://github.com/sebastian-zieba/PACMAN/blob/master/docs/source/media/Pacman_V1.gif)

# PACMAN

**Welcome to PACMAN**

ALERT: Project PACMAN is currently under heavy development. Use at your own risk.

## Installation

TODO: Add short installation guide

### With Git/GitHub

You can install PACMAN directly from source on [GitHub](https://github.com/sebastian-zieba/PACMAN) in one of two ways:
	- On the [GitHub website](https://github.com/sebastian-zieba/PACMAN), click on **Code** and **Download ZIP** followed by unpacking the distribution by opening up a terminal and typing:

		```bash
		unzip PACMAN-master.zip
		```

OR, clone the repository using ``git`` by typing:

		```bash
		git clone https://github.com/sebastian-zieba/PACMAN.git
		```

[comment]: <> (2. Navigate into the newly created directory and **install** PACMAN by running ``setup.py``:)

[comment]: <> (	```bash)

[comment]: <> (	python setup.py install)

[comment]: <> (	```)

[comment]: <> (3. Install additional **requirements** for the package by typing:)

[comment]: <> (	```bash)

[comment]: <> (	pip install -r requirements.txt)

[comment]: <> (	```)


## Documentation

Check out the docs [here](https://pacmandocs.readthedocs.io/en/latest/).




[comment]: <> (## PACMAN Steps:)



[comment]: <> (- ./wfc3_reduction/reduction/s00_table.py)

[comment]: <> (  - Creates a new run directory with the following form, eg.: ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/)

[comment]: <> (  - Copy and pastes the control file &#40;obs_par.ecf&#41; and the fit parameters file &#40;fit_par_new2.txt&#41; into the new directory)

[comment]: <> (  - Reads in all fits files and creates a table which will be saved in filelist.txt. The following information will be listed in this table: filenames filter/grism ivisit iorbit t_mjd t_visit t_orbit scan exp)

[comment]: <> (Meta data will be saved after every step in a file call something like ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/WFC3_L-98-59_Hubble15856_Meta_Save.dat)

[comment]: <> (Meta data includes information like fits headers containing observation specific parameters.  )

[comment]: <> (#)

[comment]: <> (- ./wfc3_reduction/reduction/s01_horizons.py)

[comment]: <> (  - Retrieves vector data of Hubble from JPL's HORIZONS system on https://ssd.jpl.nasa.gov/horizons_batch.cgi &#40;see Web interface on https://ssd.jpl.nasa.gov/horizons.cgi&#41;)

[comment]: <> (  - txt file with HST positions in space will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/ancil/horizons)

[comment]: <> (#)

[comment]: <> (- ./wfc3_reduction/reduction/s02_barycorr.py)

[comment]: <> (  - performs barycorr based on the t_mjd in filelist.txt. Adds another column called t_bjd)

[comment]: <> (  - Plots will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/ancil/horizons)

[comment]: <> (#)

[comment]: <> (- ./wfc3_reduction/reduction/s03_refspectra.py)

[comment]: <> (  - Downloads the bandpass of G141 or G102 and a stellar spectrum using the python package pysynphot. The product of these two will be used as a reference spectrum and for the wavelength calibration.)

[comment]: <> (  - Plots and data will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/ancil/bandpass)

[comment]: <> (#)

[comment]: <> (- ./wfc3_reduction/reduction/s10_direct_images.py)

[comment]: <> (  - The direct images will be read in and the centroid position of the star will be determined.)

[comment]: <> (  - This information will be saved in a new file called ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/xrefyref.txt)

[comment]: <> (  - A plot of the direct image will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/figs/images)

[comment]: <> (#)

[comment]: <> (- ./wfc3_reduction/reduction/s20_extract.py)

[comment]: <> (  - Creates a new directory called ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/extracted_lc/2021-08-03_12:22/)

[comment]: <> (  - The white lightcurve will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/extracted_lc/2021-08-03_12:22/lc_white.txt)

[comment]: <> (  - The spectroscopic lightcurve data will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/extracted_lc/2021-08-03_12:22/lc_spec.txt &#40;the actual binning happens in step s22&#41;)

[comment]: <> (  - Several plots will be saved in ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/figs/images)

[comment]: <> (#)

[comment]: <> (- ./wfc3_reduction/reduction/s21_ld_inputmaker.py OPTIONAL)

[comment]: <> (  - Doesnt really completely work yet! Will create limb darkening parameters if the user doesnt want to fit for them but fix them to theoretical values.)

[comment]: <> (#)

[comment]: <> (- ./wfc3_reduction/reduction/s22_bin_spectroscopic_lc.py)

[comment]: <> (  - Creates a new directory called ./run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/extracted_sp/bins12_2021-08-09_04-36-06/ and saves the binned spectroscopic lcs there. )

[comment]: <> (#)

[comment]: <> (- ./wfc3_reduction/reduction/s30_run.py)

[comment]: <> (  - Fits the spectroscopic or white lcs. Set fit_par_new2.txt for that!)

[comment]: <> (## PACMAN control file:)

[comment]: <> (This file located in ./run has to be set before running any step! Let's have a look at the different parameters here:)

[comment]: <> (Some short explanations:)

[comment]: <> (### 00)

[comment]: <> (- path     /home/zieba/Desktop/Data/L-98-59_Hubble15856/ #location of fits files)

[comment]: <> (- suffix             ima )

[comment]: <> (- which_visits  [0] # just use the 0th visit)

[comment]: <> (### 10)

[comment]: <> (- di_rmin  270 #estimated position of the star in the direct image. that cutout will be used to determine the centroid position then)

[comment]: <> (- di_rmax  290)

[comment]: <> (- di_cmin  245)

[comment]: <> (- di_cmax  265)

[comment]: <> (### 20)

[comment]: <> (- window   10: plus and minus the size of aperture in pixels of the row with the highest gradient in flux)

[comment]: <> (- diagnostics        False)

[comment]: <> (- background_thld  1000 #background threshold in counts. will just take the median below this value to determine the background flux)

[comment]: <> (- opt_extract  True # want to use optimal extraction? If not the flux in the box aperture will be simply added up)

[comment]: <> (- sig_cut  15  #optimal extraction, for cosmic rays etc)

[comment]: <> (- nsmooth  9         #optimal extraction, created smoothed spatial profile, medial smoothing filter)

[comment]: <> (- rdnoise  22.0 #constant for HST)

[comment]: <> (- output   True)

[comment]: <> (- correct_wave_shift   True)

[comment]: <> (### 22)

[comment]: <> (- wvl_min  1.125 #start of wavenlegth range to consider)

[comment]: <> (- wvl_max  1.65 #end of wavenlegth range to consider)

[comment]: <> (- wvl_bins   [6] #6 bins. user can also enter an array here like [6,12,18]. This will create three different directories with the different spec lcs then )

[comment]: <> (##30)

[comment]: <> (- toffset  2458366 #time offset which will be subtracted from the pobservations to prevent a floating point precision problem)

[comment]: <> (- run_verbose         True)

[comment]: <> (- run_output          False)

[comment]: <> (- run_show_plot       True)

[comment]: <> (- run_mcmc            True  #Do MCMC with emcee)

[comment]: <> (- run_nested          False #Do NS with dynesty)

[comment]: <> (- run_lsq             True  #Do least square &#40;with MPFit&#41; should always be true)

[comment]: <> (- run_plot_raw_data   True  #Create a plot with the raw light curve)

[comment]: <> (- run_fit_white       True  #SET THIS TRUE IF YOU WANT TO DO A WHITE LC FIT. IF YOU WANT TO FIT THE SPECTR LCS SET FALSE )

[comment]: <> (- run_divide_white    False # i usually set this false. uses the divide white method. not tested)

[comment]: <> (- run_files           ['/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-16_22-31-42_L-98-59_Hubble15856/extracted_sp/bins12_2021-08-17_23-34-29'] # example for spec fit)

[comment]: <> (- run_files           ['/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-21_02-00-29_L-98-59_Hubble15856/extracted_lc/2021-08-21_02:01/lc_white.txt'] # example for white fit)

[comment]: <> (- run_myfuncs         ['constant','upstream_downstream','model_ramp','polynomial1','transit'] #['constant', 'upstream_downstream', 'ackbar', 'polynomial1', 'transit']# ['constant', 'upstream_downstream', 'model_ramp', 'polynomial2', 'transit'])


[comment]: <> (- save_allan_plot    False)

[comment]: <> (####mcmc)

[comment]: <> (- run_nsteps          12000)

[comment]: <> (- run_nwalkers	     75)

[comment]: <> (- run_nburn           6000)

[comment]: <> (####Ns)

[comment]: <> (- run_dlogz     25000)

[comment]: <> (- run_nlive     100)


[comment]: <> (- run_clipsigma       4 #will remove outliers iteratively and fit the lc. in this case 4 sigma outliers. it will iterate a max of 4 times. i believe that works just for the white lc. might be good to set both to 0 for now)

[comment]: <> (- run_clipiters       4)


