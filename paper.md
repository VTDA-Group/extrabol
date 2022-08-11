---
title: 'ExTraBol: '
tags:
  - Python
  - astronomy
  - supernovae
  - LSST
authors:
  - name: Ian M. Thornton
    orcid: 0000-0003-2498-3008
    equal-contrib: true
    affiliation: "Penn State University" # (Multiple affiliations must be quoted)


date: 9 August 2022
bibliography: paper.bib

---

# Summary

When an observer wishes to take measurements of a supernova, or any other astronomical event, they have a couple options. The first option is using a spectrometer which spreads out the incoming photons into a spectrum of wavelengths. This method yields a lot of useful information, but it takes a significant amount of time. A quicker (and thus cheaper) option is to utilize filters to observe the event in a handful of finite wavelength windows over time. This yields broadband photometric data. The goal of `extrabol` is to extract valuable information from this broadband data in a way that is quick and user-friendly. It does this in 3 main steps:

1. Read in a data file holding observations of a supernova or other transient event over time, through any number of filters.

2. Interpolate this data using a Gaussian Process, a non-physical, statistical model utilizing covariance.

3. Fit a series of blackbody light curves to the interpolated data to estimate the bolometric luminosity, temperature, and radius of the transient over time.

`extrabol` also includes other auxiliary functions such as de-redshifting data, correcting for extinction due to galactic dust, and even converting other common data files to an extrabol-compatible format.

Because not all observations are the same, `extrabol` was built with flexibility in mind. If desired, the user can specify a known redshift, a known amount of extinction, a supernova template to assist the Gaussian Process, a minimum signal-to-noise ratio for observations, the time window to be analyzed, and more, as specified in the README. In addition to the flexibility that allows `extrabol` to work for any single object or event, it can easily be run on large batches.  


# Statement of need

Currently, about [FIND FIGURE] supernovae are observed every year. These observations must be processed by astronomers. Some tools exist already to do this such as [DISCUSS AVAILABLE OPTIONS]. In 2024, the Vera Rubin Observatory in Chile is expected to begin its Legacy Survey of Space and Time (LSST) which could increase supernovae observations by [FIND FIGURE]. Such an increase in data generates a need for new, more efficient tools to process said data. This was the motivation behind `extrabol`. When run in parallel, `extrabol` has the capability of distilling thousands of data files into useable plots in just a couple of hours [CONFIRM THIS CLAIM]. If an object of interest is found from this distillation, `extrabol` retains the flexibility to conform to unique cases.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References