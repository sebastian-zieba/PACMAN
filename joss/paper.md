---
title: '`PACMAN`: A pipeline to reduce and analyze Hubble Wide Field Camera 3 IF Grism data'
tags:
  - HST
  - python
  - astronomy
  - exoplanets
  - spectroscopy
  - photometry
authors:

  - name: Sebastian Zieba #^[zieba@mpia] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-0562-6750
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Laura Kreidberg
    orcid: 0000-0003-0514-1147
    affiliation: 1
affiliations:
 - name: Max-Planck-Institut für Astronomie, Königstuhl 17, D-69117 Heidelberg, Germany
   index: 1
 - name: Leiden Observatory, Leiden University, Niels Bohrweg 2, 2333CA Leiden, The Netherlands
   index: 2
date: xxx
bibliography: paper.bib

---

# Summary

The Hubble Space Telescope (HST) has become the preeminent workhorse facility for the characterization of extrasolar planets.
Launched in 1990 and never designed for the observations of exoplanets, the STIS spectrograph on HST was in fact used in 2002 to detect the first atmosphere ever discovered on a planet outside of our solar system [@Charbonneau2002].
With the deactivation of the Spitzer Space Telescope in 2020, HST has the two most powerful tools in space to characterize exoplanets over a broad spectral range:
(1) The Space Telescope Imaging Spectrograph (STIS; installed in 1997) in the UV and the Wide Field Camera 3 (WFC3; installed in 2009) in the Near Infrared (NIR).
With the introduction of a spatial scan mode on WFC3 [@McCullough2012] where the star is being moved perpendicular to the dispersion direction during an exposure, WFC3 observations have become very efficient due to the reduction of overhead time and the possibility of longer exposures without saturation.

For exoplanet characterization, WFC3 is used for transit and secondary eclipse spectroscopy and phase curves observations.
The instrument has two different grisms: G102 with a useful spectral range from 800 nm to up to 1150 nm and G141 encompassing 1075 nm to about 1700 nm.
The spectral range of WFC3/G141 contains several molecular species with the strongest being water (H2O) at approximately 1.4 microns leading to the successful detection of H2O in the atmosphere of over a dozen of exoplanets [@Deming2013; @Huitson2013; @Birkby2013; @Fraine2014; @Kreidberg2014b; @Evans2016].
The bluer part of WFC3, the G102 grism, has been used less often in observing programs but most notably led to the detection of Helium in the atmosphere of an exoplanet [@Spake2018].

HST will stay the most powerful space based tool for the characterization of exoplanets until the first data of the recently launched James Webb Space Telescope (JWST) reaches the machines of the observers.
Even after then, HST is expected to produce even more impactful science results due to its exquisite data.

(Maybe discuss advantage of space based data? no systematics because no atmosphere, no atmosphere so no absorption in uv and water bands and no strong IR background)

Here we present `PACMAN`, an end-to-end pipeline developed to reduce and analyze HST/WFC3 data.
The foundation of the pipeline has been already used in numerous publications [e.g.,; @Kreidberg2014a; @Kreidberg2018] and these papers have already accumulated hundreds of citations.


# Statement of need

Validating your result with an independent pipeline is always a good idea..


# Outline of the pipeline steps

The pipeline starts with the _ima_ data products provided by the Space Telescope Science Institute which can be easily accessed from [MAST](https://mast.stsci.edu/search/hst).
These files created by the WFC3 calibration pipeline, `calwf3`, have already several calibrations applied (dark subtraction, and linearity correction, flat-fielding) to each readout of the IR exposure.

In the following we want to highlight several steps in the reduction and fitting stages of the code which are typical for HST/WFC3 observations:

- **Wavelength calibration**: We create a reference spectrum out of the throughput of the respective grism (G102 or G141) and a stellar model.
The user can decide if they want to download a stellar spectrum from MAST for that or just use a black body spectrum.
This template will be then used for the wavelength calibration. 
We also determine the position of the star in the direct images which are often taken at the start of HST orbits to create an initial guess for the wavelength solution using the know dispersion of the grism.
Using the reference spectrum as a template, we determine a shift and scaling in wavelength-space that minimizes the difference between it and the first spectrum in the visit.
This first exposure in the visit in then used as the template for the following exposures in the visit.

- **Optimal extraction and outlier removal**: `PACMAN` uses an optimal extraction algorithm as presented in @Horne1986 which iteratively masks bad pixels in the image. 
We also mask bad pixels that have been flagged by `calwf3` with data quality DQ = 4 or 512\footnote{for a list of DQ flags see [https://wfc3tools.readthedocs.io/en/latest/wfc3tools/calwf3.html#data-quality-initialization-dqicorr](https://wfc3tools.readthedocs.io/en/latest/wfc3tools/calwf3.html#data-quality-initialization-dqicorr).

- **Scanning of the detector**: The majority of exoplanetary HST/WFC3 observations use the spatial scanning technique [@McCullough2012] which spreads the light perpendicular to the dispersion direction during the exposure enabled longer integration times before saturation.
The _ima_ files taken in this observation mode consist out of a number of nondestructive reads, also known as up-the-ramp samples, each of which we treat as an independent subexposure.
\autoref{fig:figure1} (left panel) shows an example for the last subexposure when using the spatial scanning together with the theoretical position of the trace.

  - **Fitting models**: `PACMAN` contains several functions to fit models which are commently used with HST data. Here some examples for the currently implemented systematic and astrophysical models.
    - systematic models:
      - visit-long polynomials
      - orbit-long exponential ramps due to charge trapping: NIR detectors like HST/WFC3 can trap photoelectrons (Smith et al. 2008) which will cause the the number of
  recorded photoelectrons that will increase exponentially creating typical hook like features
    - systematic models:
      - transit and secondary eclipse curves as implemented in `batman`
      - sinusoids for phase curve fits
      - a constant offset which accounts for the upstream-downstream effect caused by forward and reverse scanning

  The user can fit models like in autoref{eq:equation1} to the white light curve or to spectroscopic light curves. For the latter, the user can freely set the amount and locations of the bins. 
  \autoref{fig:figure1} (left panel) shows the resulting 1D spectrum and a user defined binning.
  
  \begin{equation}
  \label{eq:equation1}
  F(t) = T(t) \, (c\,S(t) + k\,t_{\rm{v}}) \, (1 - \exp(-r_1\,t_{\rm{orb}} - r_2 )),
  \end{equation}
  with T(t) being the transit model, ...

![Left panel: raw 2D spectrum. Right panel: 1D spectrum after the use of optimal extraction.\label{fig:figure1}](figures/figure1.png "title-2"){ width=99% }




![Left panel: raw spectroscopic light curve. Right panel: light curve with the best astrophysical model fit.\label{fig:figure2}](figures/figure2.png "title-2"){ width=99% }





# Dependencies

`PACMAN` uses typical dependencies of astrophysical python codes: `numpy` [@numpy2020], `matplotlib` [@matplotlib2007], `scipy` [@scipy2020] and `astropy` [@astropy2013; @astropy2018].

Other dependencies which might be required for the fitting stage depending on the model and sampler being run are: `batman` [@Kreidberg2015], `emcee` [@Foreman-Mackey2013], `dynesty` [@Speagle2020] and `corner` [@corner2016].

For the barycentric correction, `PACMAN` accesses the [API to JPL's Horizons system](https://ssd-api.jpl.nasa.gov/obsolete/horizons_batch_cgi.html).

If the user decides to use a stellar spectrum for the wavelength calibration, `PACMAN` will download the needed fits file from the [REFERENCE-ATLASES HLSP](https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/) hosted on the MAST archive.[@STScI2013]


# Documentation

The documentation for `PACMAN` can be found at [pacmandocs.readthedocs.io](https://pacmandocs.readthedocs.io/en/latest/) hosted on [ReadTheDocs](https://readthedocs.org/).
It includes most notably, a full explanation of every parameter in the _pacman control file_ (pcf), the API and an example on how to download, reduce and analyse observations of GJ 1214 b taken with HST/WFC3/G141.

# Similar tools

The only other end-to-end open source pipeline\footnote{as far as we are aware of} specifically developed for the reduction and analysis of HST/WFC3 data is [`Iraclis`](https://github.com/ucl-exoplanets/Iraclis).
For a full explanation of the different steps in this pipeline we refer to @Tsiaras2016.

Another open-source pipeline which has been used as an independent check of recent results presented in @Mugnai2021 and @Carone2021 is [`CASCADe`](https://jbouwman.gitlab.io/CASCADe/) (Calibration of trAnsit Spectroscopy using CAusal Data).
The pipeline has been applied to both Hubble and Spitzer datasets and uses causal connections within a dataset to model both transit signal and systematics.
For a more detailed discussion of `CASCADe` see Appendix 1 in @Carone2021.


# Future work

Additional fitting models are planned to be added in the future like phase curves using the SPIDERMAN package. 

(Could it do old data sets in staring mode?)

# Acknowledgements

We acknowledge B. Zawadzki for the creation of the PACMAN logo.



# References
