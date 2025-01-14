[![DOI](https://zenodo.org/badge/283353973.svg)](https://zenodo.org/badge/latestdoi/283353973)


# extrabol

`extrabol` is a Python 3.x package for rapidly and systematicaly estimating the bolometric luminosity and black body properties of (thermal) extragalactic transients from broadband UVONIR data. `extrabol` is broken down into three main steps:

- Read in and pre-process a data file holding observations of a supernova or other transient event over time, through any number of broadband filters.
- Interpolate this data using a Gaussian Process(GP), a non-physical, statistical model utilizing covariance.
- Fit a series of blackbody curves to the interpolated data to estimate the bolometric luminosity, temperature, and radius of the transient over time.

Author: Ian Thornton

## Installation and Documentation

Install with pip:

```bash
pip install extrabol
```

The latest documentation is available [here](https://extrabol.readthedocs.io/en/latest/).

## Usage

```
extrabol 'filename.dat' ARGUMENTS
```
Optional Arguments:
```
--verbose
    Increase output verbosity

-m MEAN, --mean MEAN
    Use a supernova template as the mean function for the GP; Choose \'1a\',\'1bc\', \'2l\', \'2p\', or \'0\' (for no template)
    Default = 0

-t, --show-template
    Shows supernova template on plots as dotted lines

-d DIST, --dist DIST
    Object luminosity distance in Mpc

-z REDSHIFT, --redshift REDSHIFT
    Object redshift
    If no argument is provided, redshift will be read from the input file

--dm DM
    Object distance modulus

--plot BOOL
    Output plots, Default is TRUE

--outdir OUTDIR
    A file location for outputs to be written to

--ebv EBV
    Milky Way extinction E(B-V) value, if known
    If no argument is probided the extinction will be read from the input file

--hostebv HOSTEBV
    Host extinction E(B-V)

-s START, --start START
    Start time of analysis window relative to peak luminosity
    Default = -50 days

-e END, --end END
    End time of analysis window relative to peak luminosity
    Default = 150 days

-snr SNR
    Minimum signal to noise ratio for observations
    Default = 4.0

-wc, --wvcorr
    Use the redshift-corrected wavelength values for extinction calculations

-mc, --use-mcmc
    Use a Markov Chain Monte Carlo to fit black bodies instead of curve_fit.
    This provides better error estimates but takes much longer.

--T_max T_MAX
    Modify the prior on temperature for blackbody fits by specifying a maximum temperature.
    Default = 40,000K
```
## Input Files

Inputs to extrabol must be .dat files that conform to the following format.

The first two lines must contain redshift and Milky Way extinction E(b-v) respectively. If these values are unknown, simply put 0.0.
The following lines contain observational data in 5 columns as shown below:

```
Time(MJD)   Apparent Magnitude   Error(in magnitudes)   Filter SVO ID   Type of magnitude (AB or Vega)
```
Any white space can be used as the column delimiter. NaNs, non-detections, and data points with no error bars should not be included.
An example input file can be found under extrabol/example.

## Example Input

```python
extrabol ./example/PSc000174_extrabol.dat --verbose -m 1a
```

This example may take a minute the first time you run it as astroquery fetches filter information!

Filters must be specified by their [SVO](http://svo2.cab.inta-csic.es/svo/theory/fps3/) ID.
