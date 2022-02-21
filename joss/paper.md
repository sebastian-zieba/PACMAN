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

Here we present PACMAN, a pipeline capable of reducing and analysing HST/WFC3 data. 
The foundation of the pipeline has been already used in numerous publications [e.g.,; @Kreidberg2014a; @Kreidberg2018] and these papers have already accumulated hundreds of citations.


# Statement of need

Hubble data reduction aint ez clap

In the following we want to highlight several steps in the reduction and fitting stages of the code which are typical for HST/WFC3 observations:

Reduction:

- Wavelength calibration: stellar and bandpass template

- Outlier removal (caused for example by cosmic rays) is done by using flags already existent in the ima fits files and optimal extraction [@Horne1986].

- scanning [@McCullough2012]

Fitting:

- Charge trapping: Fit exponential [@Zhou2017] NIR detectors like HST/WFC3 can trap photoelectrons (Smith et al. 2008) which will cause the the number of
recorded photoelectrons that will increase exponentially creating typical hook like features.

- scanning [@McCullough2012]: constant offset

- divide white vs xxx (the other one)




\autoref{fig:example}

![Left panel: raw 2D spectrum. Right panel: 1D spectrum after the use of optimal extraction.\label{fig:figure1}](figures/figure1.png "title-2"){ width=99% }

![Left panel: raw spectroscopic light curve. Right panel: light curve with the best astrophysical model fit.\label{fig:figure2}](figures/figure2.png "title-2"){ width=99% }





# Dependencies

PACMAN uses typical dependencies of astrophysical python codes: `numpy` [@numpy2020], `matplotlib` [@matplotlib2007], `scipy` [@scipy2020] and `astropy` [@astropy2013; @astropy2018].

Other dependencies which might be required for the fitting stage depening on the model and sampler being run are: `batman` [@Kreidberg2015], `emcee` [@Foreman-Mackey2013], `dynesty` [@Speagle2020] and `corner` [@corner2016].

For the barycentric correction, PACMAN will access the [API to JPL's Horizons system](https://ssd-api.jpl.nasa.gov/obsolete/horizons_batch_cgi.html).

If the user decides to use a stellar spectrum for the wavelength calibration, PACMAN will download the needed fits file from the [REFERENCE-ATLASES HLSP](https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/) hosted on the MAST archive.[@STScI2013]


# Documentation

The documentation for PACMAN can be found [here](https://pacmandocs.readthedocs.io/en/latest/) hosted on Read The Docs.
It includes most notably, a full explanation of every parameter in the pacman control file (pcf), the API and an example on how to download, reduce and analyse observations of GJ 1214 b taken with HST/WFC3/G141.

# Similar tools

`Iraclis` [@Tsiaras2016]

# Future work



# Acknowledgements

We acknowledge contributions from ....

# References
