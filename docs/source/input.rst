Input Format
============

Inputs to extrabol must be .dat files that conform to the following format.

The first two lines must contain redshift and Milky Way extinction E(b-v) respectively. If these values are unknown, simply put 0.0. The following lines contain observational data in 5 columns:

1. Time(MJD)
2. Apparent Magnitude
3. Error(in magnitudes)
4. Filter SVO ID
5. Type of magnitude (AB or Vega)

Any white space can be used as the column delimiter. NaNs, non-detections, and data points with no error bars should not be included. An example input file can be found under ``extrabol/example``.