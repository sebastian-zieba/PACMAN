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

pipeline explained: @Kreidberg2014 @Kreidberg2018 

# Statement of need

complicated... 
scanning @McCullough2012
optimal extraction @Horne1986
charge trapping @Zhou2017

# Dependencies

`numpy` [@numpy2020]
`matplotlib` [@matplotlib2007]
`scipy` [@scipy2020]
`astropy` [@astropy2013; @astropy2018]


optional:
`batman` [@Kreidberg2015]
`emcee` [@Foreman-Mackey2013]
`dynesty` [@Speagle2020]
`corner` [@corner2016]

# Documentation

# Similar tools

`Iraclis` [@Tsiaras2016]

# Future work

# How to cite




# Mathematics

$f(x) = e^{\pi/x}$

\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.


# Figures

Figure sizes can be customized by adding an optional second parameter: \autoref{fig:example} h

![alt-text-1](figures/di_0.png "title-1"){ width=30% } ![alt-text-2](figures/trace_61.png "title-2"){ width=40% }

![alt-text-3](figures/spec_bins12.png "title-1"){ width=30% } ![alt-text-4](figures/raw_lc_0.png "title-2"){ width=30% } ![alt-text-5](figures/fit_lc_0_2022-02-15_22-34-53.png "title-2"){ width=30% }

![Caption for example figure.\label{fig:example}](figures/joss-logo.png)


# Acknowledgements

We acknowledge contributions from ....

# References
