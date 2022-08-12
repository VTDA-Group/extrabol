Get bolometric estimates of SNe

Note: You need the absolute LATEST version of astroquery to work with SVO

Author: Ian Thornton

Install with pip:

```bash
pip install extrabol
```
# Usage

```
extrabol 'filename.dat' ARGUMENTS
```
Optional Arguments:
```
--verbose
    Increase output verbosity

-m MEAN, --mean MEAN
    Use a supernova template as the mean function for the GP; Choose \'1a\',\'1bc\', \'2l\', \'2p\', or \'0\' (for no template)

-t, --show_template
    Shows supernova template on plots as dotted lines

-d DIST, --dist DIST
    Object luminosity distance

-z REDSHIFT, --redshift REDSHIFT
    Object redshift

--dm DM
    Object distance modulus

--plot BOOL
    Output plots, Default is TRUE

--outdir OUTDIR
    A file location for outputs to be written to

--ebv EBV
    MWebv value, if known

--hostebv HOSTEBV
    Host B-V

-s START, --start START
    Start time of analysis window

-e END, --end END
    End time of analysis window

-snr SNR
    Minimum signal to noise ratio for observations

-wc, --wvcorr
    Use the redshift-corrected wavelength values for extinction calculations

--use_mcmc
    Use a Markov Chain Monte Carlo to fit black bodies instead of curve_fit.
    This provides better error estimates but takes much longer.
```
# Input Files

Inputs to extrabol must be .dat files that conform to the following format.

The first two lines must contain redshift and MWebv respectively. If these values are unknown, simply put 0.0.
The following lines contain observational data in 5 columns as shown below:

```
Time(MJD)   Apparent Magnitude   Error(in magnitudes)   Filter SVO ID   Type of magnitude (AB or Vega)
```
An example input file can be found under extrabol/example

# Example input

```python
extrabol ./example/PSc000174_extrabol.dat --verbose -m 1a
```

This example may take a minute the first time you run it as astroquery fetches filter information!

Filters must be specified by their [SVO](http://svo2.cab.inta-csic.es/svo/theory/fps3/) ID.

