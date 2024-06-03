Usage
============

extrabol can be used directly from commandline:

.. code-block:: bash

    extrabol 'filename.dat' ARGUMENTS


The following arguments are available: 

.. code-block:: bash

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

    -mc, --use-mcmc
        Use a Markov Chain Monte Carlo to fit black bodies instead of curve_fit.
        This provides better error estimates but takes much longer.

    --T_max T_MAX
        Modify the prior on temperature for blackbody fits by specifying a maximum temperature.
        Default = 40,000K

    -k, --kernel-width
        The width (:math:`r^2`) of the GP kernel in the (time, wavelength) direction.
        If not given, the kernel width will be optimized.
