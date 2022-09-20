Get bolometric estimates of SNe

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

-t, --show-template
    Shows supernova template on plots as dotted lines

-d DIST, --dist DIST
    Object luminosity distance in Mpc

-z REDSHIFT, --redshift REDSHIFT
    Object redshift

--dm DM
    Object distance modulus

--plot BOOL
    Output plots, Default is TRUE

--outdir OUTDIR
    A file location for outputs to be written to

--ebv EBV
    Milky Way extinction E(B-V) value, if known

--hostebv HOSTEBV
    Host extinction E(B-V)

-s START, --start START
    Start time of analysis window relative to peak luminosity

-e END, --end END
    End time of analysis window relative to peak luminosity

-snr SNR
    Minimum signal to noise ratio for observations

-wc, --wvcorr
    Use the redshift-corrected wavelength values for extinction calculations

--use-mcmc
    Use a Markov Chain Monte Carlo to fit black bodies instead of curve_fit.
    This provides better error estimates but takes much longer.
```
# Input Files

Inputs to extrabol must be .dat files that conform to the following format.

The first two lines must contain redshift and Milky Way extinction E(b-v) respectively. If these values are unknown, simply put 0.0.
The following lines contain observational data in 5 columns as shown below:

```
Time(MJD)   Apparent Magnitude   Error(in magnitudes)   Filter SVO ID   Type of magnitude (AB or Vega)
```
Any white space can be used as the column delimiter. NaNs, non-detections, and data points with no error bars should not be included.
An example input file can be found under extrabol/example.

# Example input

```python
extrabol ./example/PSc000174_extrabol.dat --verbose -m 1a
```

This example may take a minute the first time you run it as astroquery fetches filter information!

Filters must be specified by their [SVO](http://svo2.cab.inta-csic.es/svo/theory/fps3/) ID.

