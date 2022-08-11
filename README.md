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
```
# Example input

```python
extrabol ./example/PSc000174_extrabol.dat --verbose -m 1a
```

This example may take a minute the first time you run it as astroquery fetches filter information!

Filters must be specified by their [SVO](http://svo2.cab.inta-csic.es/svo/theory/fps3/) ID.

