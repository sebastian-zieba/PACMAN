.. _Outputs:

Outputs
=======

Here we summarize all the files created by PACMAN when doing a complete run from s00 to s30.


Stage 00
''''''''

directories
-----------
- workdir (e.g.: run_2022-12-16_10-46-31_dec2022) 
    the main directory of your analysis

- workdir/figs
    all figures will be saved into here

- workdir/figs/s00_obs_dates
    

plots
------------
- workdir/figs/s00_obs_dates/obs_dates.png
- workdir/figs/s00_obs_dates/obs_dates_all.png

other files
------------

- workdir/WFC3_dec2022_Meta_Save.dat
    the meta file: this file will collect a lot of useful information when finishing a stage and then accessed by the next stage.

- workdir/filelist.txt
    a list containing all *ima.fits files in the datadir and some information from their headers.

- workdir/fit_par.txt
    this fit_par.txt will be used for the fitting in s30

- workdir/obs_par.pcf
    this is the main obs_par.pcf file for all future stages

Stage 01
''''''''

directories
-----------
- ancil
- ancil/horizons

other files
------------

- ancil/horizons/horizons_results_v0.txt


Stage 02
''''''''

directories
-----------
- workdir/figs/s02_barycorr

plots
------------
- workdir/figs/s02_barycorr/bjdcorr_horizons_results_v0.png


Stage 03
''''''''

directories
-----------

- workdir/ancil/refspec

- workdir/ancil/stellar_models

- workdir/figs/s03_refspec

- workdir/figs/s03_smooth

plots
------------

- workdir/figs/s03_refspec/refspec.png

- workdir/figs/s03_smooth/smooth.png

other files
------------

- workdir/ancil/stellar_models/k93models/kp03_3500.fits



Stage 10
''''''''

directories
-----------

- workdir/figs/s10_images

plots
------------

- workdir/figs/di_0.png
- workdir/figs/quick_di0.png

other files
------------

- workdir/xrefyref.txt




Stage 21
''''''''

directories
-----------

- workdir/extracted_lc/2022-12-16_11-22-12
- workdir/figs/s20_badmask
- workdir/figs/s20_bkg_evo
- workdir/figs/s20_bkg_hist
- workdir/figs/s20_drift
- workdir/figs/s20_optextr
- workdir/figs/s20_refspec_fit
- workdir/figs/s20_sp1d
- workdir/figs/s20_sp1d_diff
- workdir/figs/s20_sp2d
- workdir/figs/s20_trace
- workdir/figs/s20_utr
- workdir/figs/s20_utr_aper_evo

plots
------------

- workdir/figs/s20_badmask/badmask_0
- workdir/figs/s20_bkg_evo/bkg_evo.png
- workdir/figs/s20_bkg_hist/bkg_hist0-0.png
- workdir/figs/s20_drift/drift.png
- workdir/figs/s20_optextr/optextr0-0.png
- workdir/figs/s20_refspec_fit/refspec_fit_0.png
- workdir/figs/s20_sp1d/sp1d_0.png
- workdir/figs/s20_sp1d_diff/sp1d_diff_0.png
- workdir/figs/s20_sp2d/sp2d_0.png
- workdir/figs/s20_trace/trace_0.png
- workdir/figs/s20_utr/utr0-0.png
- workdir/figs/s20_utr_aper_evo/utr_aper_evo.png

other files
------------

- workdir/extracted_lc/2022-12-16_11-22-12/diagnostics.txt
- workdir/extracted_lc/2022-12-16_11-22-12/lc_spec.txt
- workdir/extracted_lc/2022-12-16_11-22-12/lc_white.txt
- workdir/extracted_lc/2022-12-16_11-22-12/obs_par.pcf
    this pcf is just a copy of the pcf which was used when s20 was run. The purpose of this copy is just so that the user can trace back which settings were uses when running this stage.


Stage 20
''''''''

directories
-----------

- workdir/extracted_sp
- workdir/extracted_sp/bins11_2022-12-16_11-52-46

plots
------------

- workdir/extracted_sp/bins11_2022-12-16_11-52-46/spec_bins11.png

other files
------------

- workdir/extracted_sp/bins11_2022-12-16_11-52-46/speclc1.158.txt
- workdir/extracted_sp/bins11_2022-12-16_11-52-46/wvl_table.dat


Stage 30 (white)
''''''''''''''''

directories
-----------

- workdir/fit_white
- workdir/fit_white/fit_2022-12-16_12-02-03_dec2022
- workdir/fit_white/fit_2022-12-16_12-02-03_dec2022/fit_lc
- workdir/fit_white/fit_2022-12-16_12-02-03_dec2022/lsq_res
- workdir/fit_white/fit_2022-12-16_12-02-03_dec2022/mcmc_res
- workdir/fit_white/fit_2022-12-16_12-02-03_dec2022/raw_lc

plots
------------

- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/fit_lc/lsq_lc_bin0_wvl1.400.png
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/fit_lc/mcmc_lc_bin0_wvl1.400.png
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/fit_lc/newfit_lc_bin0_wvl1.400.png

- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/lsq_res/corr_plot_bin0_wvl1.400.png
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/lsq_res/lsq_params_vs_wvl.png

- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/mcmc_res/corr_plot_bin0_wvl1.400.png
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/mcmc_res/mcmc_chains_bin0_wvl1.400.png
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/mcmc_res/mcmc_chains_noburn_bin0_wvl1.400.png
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/mcmc_res/mcmc_pairs_bin0_wvl1.400.png
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/mcmc_res/mcmc_params_vs_wvl.png

- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/raw_lc/raw_lc_bin0_wvl1.400.png

other files
------------
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/fit_lc/fit_lc_data_model_bin0_wvl1.400.txt
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/fit_lc/fit_lc_data_nosys_bin0_wvl1.400.txt

- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/lsq_res/corr_data_bin0_wvl1.400.txt
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/lsq_res/lsq_res_bin0_wvl1.400.txt

- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/mcmc_res/corr_data_bin0_wvl1.400.txt
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/mcmc_res/mcmc_out_bin0_wvl1.400.p
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/mcmc_res/mcmc_res_bin0_wvl1.400.txt
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/raw_lc/raw_lc_data_bin0_wvl1.400.txt



- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/fit_par.txt
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/fit_results.txt
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/obs_par.pcf
- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/white_systematics.txt
    file for divide white with the best systematics model after the least squares fit

- workdir/fit_white/fit_2022-12-16_14-51-31_dec2022/white_systematics_mcmc.txt
    file for divide white with the best systematics model after the mcmc


Stage 30 (spec)
'''''''''''''''

directories
-----------

- workdir/fit_spec
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/fit_lc
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/lsq_res
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/raw_lc

plots
------------


- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/fit_lc/lsq_lc_bin0_wvl1.158.png
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/fit_lc/newfit_lc_bin0_wvl1.158.png



- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/lsq_res/corr_plot_bin0_wvl1.158.png
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/lsq_res/lsq_params_vs_wvl.png

- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/lsq_res/lsq_rprs.png



- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res/corr_plot_bin0_wvl1.158.png
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res/mcmc_chains_bin0_wvl1.158.png
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res/mcmc_chains_noburn_bin0_wvl1.158.png
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res/mcmc_pairs_bin0_wvl1.158.png
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res/mcmc_params_vs_wvl.png
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res/mcmc_rprs.png


- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/raw_lc/raw_lc_bin0_wvl1.158.png

other files
------------
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/fit_par.txt
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/fit_results.txt

	Here an example of a fit_resutlts.txt file:

	.. include:: media/s30/spectroscopic/fit_results.txt

	It includes the following columns:
		- wave: mid wavelength of the bin
		- rms_pred: expected photon noise 
		- rms: rms of residuals
		- chi2: chi squared: ``= np.sum((self.resid/data.err)**2)``
		- chi2red: ``= chi2/dof``, with dof being the degrees of freedom 
		- bic: baysian information criterium: ``= -2. * ln_like + nfree_param * np.log(npoints)``
		- ln_like: log likelihood ``= (-0.5 * (np.sum((residuals/errors)**2 + np.log(2.0 * np.pi) + 2 * np.log(errors)))``
		- npoints: number of data points
		- nfree_param: number of free parameters
		- dof: degree of freedom

    The next parameters will just appear if the errorbars were rescaled in the analysis. So if either "meta.uncert_rescale = True" or "uncmulti" was used as a model. The "notrescaled" refers to the fact that for the calculation of these chi2, chi2red, bic and ln_like values, the original errorbars were NOT rescaled. 

		- chi2_notrescaled
		- chi2red_notrescaled
		- bic_notrescaled
		- ln_like_notrescaled


- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/obs_par.pcf

- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/fit_lc/fit_lc_data_model_bin0_wvl1.158.txt
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/fit_lc/fit_lc_data_nosys_bin0_wvl1.158.txt

- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/lsq_res/corr_data_bin0_wvl1.158.txt
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/lsq_res/lsq_res_bin0_wvl1.158.txt
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/lsq_res/lsq_rprs.txt

- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res/corr_data_bin0_wvl1.158.txt
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res/mcmc_out_bin0_wvl1.158.p
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res/mcmc_res_bin0_wvl1.158.txt
- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/mcmc_res/mcmc_rprs.txt

- workdir/fit_spec/fit_2022-12-16_16-50-21_dec2022/raw_lc/raw_lc_data_bin1_wvl1.204.txt









