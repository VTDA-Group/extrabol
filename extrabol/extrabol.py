#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import george
from scipy.optimize import minimize, curve_fit
import argparse
from astropy.cosmology import Planck13 as cosmo
from astropy import units as u
import os
from astropy.table import Table
from astropy.io import ascii
import matplotlib.cm as cm
import sys
from scipy import interpolate as interp
from george.modeling import Model
import extinction
import emcee
import importlib_resources
import pandas as pd
import re


# Define a few important constants that will be used later
epsilon = 0.001
c = 2.99792458E10  # cm / s
sigsb = 5.6704e-5  # erg / cm^2 / s / K^4
h = 6.62607E-27
ang_to_cm = 1e-8
k_B = 1.38064852E-16  # cm^2 * g / s^2 / K


def bbody(lam, T, R):
    '''
    Calculate BB L_lam (adapted from superbol, Nicholl, M. 2018, RNAAS)

    Parameters
    ----------
    lam : float
        Reference wavelengths in Angstroms
    T : float
        Temperature in Kelvin
    R : float
        Radius in cm

    Output
    ------
    L_lam in erg/s/cm
    '''

    lam_cm = lam * ang_to_cm
    exponential = (h*c) / (lam_cm*k_B*T)
    blam = ((2.*np.pi*h*c**2) / (lam_cm**5)) / (np.exp(exponential)-1.)
    area = 4. * np.pi * R**2
    lum = blam * area

    return lum


import pandas as pd
import re

def read_snana_file(file_path):
    '''
    Reads a SNANA file and extracts metadata and observational data.
    Written by ChatGPT! Thanks, ChatGPT.

    Parameters
    ----------
    file_path : str
        The path to the SNANA file.

    Returns
    -------
    metadata : dict
        A dictionary containing the metadata from the file. Metadata includes all lines before the "NOBS" line that lack a hashtag/pound sign.
        If a value contains a '+-' (e.g., '0.0076 +- 0.0003'), only the float value before '+-' is extracted.
    data_df : pandas.DataFrame
        A pandas DataFrame containing the observational data. The columns are defined by the 'VARLIST' line, and the observations follow the 'OBS' lines.
        Each row corresponds to an observation, with columns as specified in 'VARLIST'.
    
    Example
    -------
    >>> metadata, data_df = read_snana_file('path_to_your_snana_file.txt')
    >>> print(metadata)
    {'SURVEY': 'Archival', 'SNID': 'SN2010bc', 'IAUC': 'SN2010bc', 'RA': '162.02942', 'DECL': '56.85011', 'MWEBV': '0.0076', 'REDSHIFT_FINAL': '0.2440', 'SEARCH_PEAKMJD': '55222.5', 'FILTERS': 'griz'}
    >>> print(data_df.head())
          MJD                  FLT        MAG     MAGERR MAGTYPE
    0  55173.6  PAN-STARRS/PS1.g  26.213299   3.420972      AB
    1  55191.6  PAN-STARRS/PS1.g  29.423878  60.670588      AB
    '''

    metadata = {}
    data_lines = []
    varlist = []
    header_found = False

    with open(file_path, 'r') as file:
        for line in file:
            stripped_line = line.strip()

            if stripped_line.startswith('#'):
                continue

            if stripped_line.startswith('NOBS'):
                header_found = True
                continue

            if not header_found:
                if ':' in line:
                    key, value = line.split(':', 1)
                    value = value.strip()
                    # Check for +- and extract only the float value before the +-
                    if '+-' in value:
                        value = value.split('+-')[0].strip()
                    # Extract float value if present
                    float_value = re.findall(r"[-+]?\d*\.\d+|\d+", value)
                    if float_value:
                        value = float_value[0]
                    metadata[key.strip()] = value
            else:
                if line.startswith('VARLIST:'):
                    varlist = line.split()[1:]
                elif line.startswith('OBS:'):
                    data_lines.append(line.split()[1:])

    # Create pandas DataFrame from observation data
    if varlist and data_lines:
        data_df = pd.DataFrame(data_lines, columns=varlist)
        # Identify columns that should remain as strings
        str_columns = ['FLT', 'MAGTYPE']
        # Convert numeric columns to appropriate data types
        for col in varlist:
            if col not in str_columns:
                data_df[col] = pd.to_numeric(data_df[col], errors='coerce')

    return metadata, data_df


def read_in_photometry(filename, dm, redshift, start, end, snr, mwebv,
                       hostebv, verbose):
    '''
    Read in SN file

    Parameters
    ----------
    filename : string
        Input file name. Filetype will be assumed from filename
    dm : float
        Distance modulus
    redshift : float
        redshift
    start : float
        The first time point to be accepted
    end : float
        The last time point to be accepted
    snr : float
        The lowest signal to noise ratio to be accepted
    mwebv : float
        Milky Way extinction to be removed (in the observed frame)
    hostebv : bool
        Host-galaxy extinction to be removed (in the rest frame)

    Output
    ------
    lc : numpy.array
        Light curve array
    wv_corr : float
        Mean of wavelengths, used in GP pre-processing
    flux_corr : float
        Flux correction, used in GP pre-processing
    my_filters : list
        List of filter names
    '''

    #if '.json' in filename:
    phases = []
    errs = []
    observed_filters = []
    observed_mags = []

    if '.snana' in filename:
        print('Assuming SNANA format')
        metadata, data_df = read_snana_file(filename)

        phases = np.array(data_df['MJD'].values, dtype=float)
        errs = np.array(data_df['MAGERR'].values, dtype=float)
        observed_filters = np.array(data_df['FLT'].values, dtype=str)
        observed_mags = np.asarray(data_df['MAG'].values, dtype=float)
        mag_types = np.asarray(data_df['MAGTYPE'].values, dtype=str)

    else:

        photometry_data = np.loadtxt(filename, dtype=str, skiprows=2)
        # Extract key information into seperate arrays
        phases = np.asarray(photometry_data[:, 0], dtype=float)
        errs = np.asarray(photometry_data[:, 2], dtype=float)
        observed_filters = np.asarray(photometry_data[:,3], dtype=str)
        observed_mags = np.asarray(photometry_data[:,1], dtype=float)
        mag_types = np.asarray(photometry_data[:,-1], dtype=str)

    filter_data = importlib_resources.files('extrabol.filter_data') / 'fps.xml'
    index = Table.read(filter_data)
    filterIDs = np.asarray(index['filterID'].data, dtype=str)
    wavelengthEffs = np.asarray(index['WavelengthEff'].data, dtype=float)
    widthEffs = np.asarray(index['WidthEff'].data, dtype=float)
    zpts_all = np.asarray(index['ZeroPoint'].data, dtype=str)

    # Extract filter names and effective wavelengths
    wv_effs = []
    width_effs = []
    my_filters = []
    for ufilt in observed_filters:
        gind = np.where(filterIDs == ufilt)[0]
        if len(gind) == 0:
            sys.exit('Cannot find ' + str(ufilt) + ' in SVO.')
        wv_effs.append(wavelengthEffs[gind][0])
        width_effs.append(widthEffs[gind][0])
        my_filters.append(ufilt)
    wv_effs = np.asarray(wv_effs)

    # Convert brightness data to flux
    zpts = []
    fluxes = []
    for i, datapoint in enumerate(observed_mags):
        mag = float(datapoint) - dm + 2.5 * np.log10(1. + redshift)  # cosmological k-correction
        if mag_types[i] == 'AB':
            zp = 0.
        else:
            gind = np.where(filterIDs == observed_filters[i])
            zp = 2.5 * np.log10(float(zpts_all[gind[0]][0]) / 3631.)
        zpts.append(zp)

        # 'fluxes' is the negative absolute AB magnitude
        # This is easier on the Gaussian Process
        flux = zp - mag
        fluxes.append(flux)

    # Remove extinction
    ext = extinction.fm07(wv_effs, mwebv * 3.1) + extinction.fm07(wv_effs / (1. + redshift), hostebv * 3.1)
    for i in np.arange(len(fluxes)):
        fluxes[i] = fluxes[i] + ext[i]

    # The GP prefers when values are relatively close to zero
    # so we adjust wavelength and flux accordingly
    # This will be undone after interpolation
    wv_corr = np.mean(wv_effs / (1.+redshift))
    flux_corr = np.min(fluxes) - 1.0
    wv_effs = (wv_effs / (1.+redshift)) - wv_corr
    fluxes = np.asarray(fluxes) - flux_corr

    # Eliminate any data points bellow threshold snr
    gis = []
    for i in np.arange(len(phases)):
        if (1/errs[i]) >= snr:
            gis.append(i)
    gis = np.asarray(gis, dtype=int)


    phases = phases[gis]
    fluxes = fluxes[gis]
    wv_effs = wv_effs[gis]
    errs = errs[gis]
    width_effs = np.asarray(width_effs)
    width_effs = width_effs[gis]
    my_filters = np.asarray(my_filters)
    my_filters = my_filters[gis]

    # Set the peak flux to t=0
    peak_i = np.argmax(fluxes)
    if verbose:
        print('Peak Luminosity occurs at MJD',phases[peak_i])
    phases = (np.asarray(phases) - phases[peak_i]) / (1 + redshift)

    # Eliminate any data points outside of specified range
    # With respect to first data point
    gis = []
    for i in np.arange(len(phases)):
        if phases[i] <= end and phases[i] >= start:
            gis.append(i)
    gis = np.asarray(gis, dtype=int)

    phases = phases[gis]
    fluxes = fluxes[gis]
    wv_effs = wv_effs[gis]
    errs = errs[gis]
    width_effs = width_effs[gis]
    my_filters = my_filters[gis]

    lc = np.vstack((phases, fluxes, wv_effs / 1000., errs, width_effs))

    return lc, wv_corr, flux_corr, my_filters


def chi_square(dat, model, uncertainty):
    '''
    Calculate the chi squared of a model given a set of data

    Parameters
    ----------
    dat : numpy.array
        Experimental data for the model to be tested against
    model : numpy.array
        Model data being tested
    uncertainty : numpy.array
        Error on experimental data

    Output
    ------
    chi2 : float
        the chi sqaured value of the model
    '''

    chi2 = 0.
    for i in np.arange(len(dat)):
        chi2 += ((model[i]-dat[i]) / uncertainty[i])**2.

    return chi2

def generate_template(filter_wv, sn_type):
    '''
    Prepare and interpolate SN1a Template

    Parameters
    ----------
    fiter_wv : numpy.array
        effective wavelength of filters in Angstroms
    sn_type : string
        The type of supernova template to be used for GP mean function

    Output
    ------
    temp_interped : RectBivariateSpline object
        interpolated template
    '''



    my_template_file = importlib_resources.files('extrabol.template_bank') / ('smoothed_sn' + sn_type + '.npz')
    template = np.load(my_template_file)
    temp_times = template['time']
    temp_wavelength = template['wavelength']
    temp_f_lambda = template['f_lambda']

    # The template is too large, so we thin it out
    # First chop off unnecessary ends
    gis = []
    for i in np.arange(len(temp_wavelength)):
        if temp_wavelength[i] < np.amax(filter_wv) and\
                temp_wavelength[i] > np.amin(filter_wv):
            gis.append(i)
    temp_times = temp_times[gis]
    temp_wavelength = temp_wavelength[gis]
    temp_f_lambda = temp_f_lambda[gis]

    # Remove every other time point
    gis = []
    for i in np.arange(len(temp_times)):
        if temp_times[i] % 2. == 0:
            gis.append(i)
    temp_times = temp_times[gis]
    temp_wavelength = temp_wavelength[gis]
    temp_f_lambda = temp_f_lambda[gis]

    # Remove every other wavelength point
    gis = []
    for i in np.arange(len(temp_wavelength)):
        if temp_wavelength[i] % 20. == 0:
            gis.append(i)
    temp_times = temp_times[gis]
    temp_wavelength = temp_wavelength[gis]
    temp_f_lambda = temp_f_lambda[gis]

    # Remove initial rise
    # If the data point is very dim, it likely has a low snr
    gis = []
    for i in np.arange(len(temp_times)):
        if temp_times[i] >= 1.:
            gis.append(i)
    temp_times = temp_times[gis]
    temp_wavelength = temp_wavelength[gis]
    temp_f_lambda = temp_f_lambda[gis]

    # Set peak flux to t=0
    peak_i = np.argmax(temp_f_lambda)
    temp_times = np.asarray(temp_times) - temp_times[peak_i]

    # RectBivariateSpline requires
    # that x and y are 1-d arrays, strictly ascending
    temp_times_u = np.unique(temp_times)
    temp_wavelength_u = np.unique(temp_wavelength)
    temp_f_lambda_u = np.zeros((len(temp_times_u), len(temp_wavelength_u)))
    for i in np.arange(len(temp_times_u)):
        gis = np.where(temp_times == temp_times_u[i])
        temp_f_lambda_u[i, :] = temp_f_lambda[gis]
    # Template needs to be converted to log(flux) to match data
    for i in np.arange(len(temp_wavelength_u)):
        wv = temp_wavelength_u[i]
        temp_f_lambda_u[:, i] = 2.5 * np.log10((wv**2) * temp_f_lambda_u[:, i])

    temp_interped = interp.RectBivariateSpline(temp_times_u, temp_wavelength_u,
                                               temp_f_lambda_u)

    return temp_interped


def fit_template(wv, template_to_fit, filts, wv_corr, flux, time,
                 errs, z, output_chi=False, output_params=True):
    '''
    Get parameters to roughly fit template to data

    Parameters
    ----------
    wv : numpy.array
        wavelenght of filters in angstroms
    template_to_fit : RectBivariateSpline object
        interpolated template
    filts : numpy.array
        normalized wavelength values for each obseration
    wv_corr : float
        Mean of wavelengths, used in GP pre-processing
    flux : numpy.array
        flux data from observations
    time : numpy.array
        time data from observations
    errs : numpy.array
        errors on flux data
    z : float
        redshift
    output_chi : bool
        If true, function returns chi squared value
    output_params : bool
        If true, function returns optimal parameters

    Output
    ------
    A_opt : float
        multiplicative constant to be applied to template flux values
    t_c_opt : float
        additive constant to line up template and data in time
    t_s_opt : float
        multiplicative constant to scale the template in the time dimention
    chi2 : float
        The chi squared value for the given parameters
    '''

    A_opt = []
    t_c_opt = []
    t_s_opt = []
    chi2 = []

    # Fit the template to the data for each filter used
    # and choose the set of parameters with the lowest total chi2
    for wavelength in wv:
        # A callable function to test chi2 later on
        def model(time, filt, A, t_c, t_s):
            time_sorted = sorted(time)
            time_corr = np.asarray(time_sorted) * 1./t_s + t_c
            mag = template_to_fit(time_corr, filt) + A  # log(flux), not mag
            mag = np.ndarray.flatten(mag)
            return mag

        # curve_fit won't know what to do with the filt param
        # so I need to modify it slightly
        def curve_to_fit(time, A, t_c, t_s):
            mag = model(time, wavelength, A, t_c, t_s)
            return mag

        # Collect the data points coresponding to the current wavelength
        gis = np.where(filts*1000 + wv_corr == wavelength)
        dat_fluxes = flux[gis]
        dat_times = time[gis]
        dat_errs = errs[gis]
        popt, pcov = curve_fit(curve_to_fit, dat_times, dat_fluxes,
                               p0=[20, 0, 1+z], maxfev=8000,
                               bounds=([-np.inf, -np.inf, 0], np.inf))
        A_opt.append(popt[0])
        t_c_opt.append(popt[1])
        t_s_opt.append(popt[2])

        # Test chi2 for this set of parameters over all filters
        param_chi = 0
        for filt in wv:
            m = model(dat_times, filt, popt[0], popt[1], popt[2])
            param_chi += chi_square(dat_fluxes, m, dat_errs)
        chi2.append(param_chi)

    # Choose the template with the minimum chi2
    gi = np.argmin(chi2)
    chi2 = chi2[gi]
    A_opt = A_opt[gi]
    t_c_opt = t_c_opt[gi]
    t_s_opt = t_s_opt[gi]

    if output_chi:
        if output_params:
            return A_opt, t_c_opt, t_s_opt, chi2
        else:
            return chi2
    else:
        if not output_params:
            return 0
        else:
            return A_opt, t_c_opt, t_s_opt


def test(lc, wv_corr, z):
    '''
    Test every available template for the lowest possible chi^2

    Parameters
    ----------
    lc : numpy.array
        LC array
    wv_corr : float
        mean of wavelength, needed to find wavelength in angstroms
    z : float
        redshift

    Output
    ------
    best_temp : string
        The name of the template with the lowest chi squared value
    '''
    lc = lc.T

    # Extract necissary information from lc
    time = lc[:, 0]
    flux = lc[:, 1]
    filts = lc[:, 2]
    errs = lc[:, 3]
    ufilts = np.unique(lc[:, 2])
    ufilts_in_angstrom = ufilts*1000.0 + wv_corr

    # Generate a template for each available supernova type
    # Then fit and test each one for lowest possible chi2
    templates = ['1a', '1bc', '2p', '2l']
    chi2 = []
    for template in templates:
        template_to_fit = generate_template(ufilts_in_angstrom, template)
        chi2.append(fit_template(ufilts_in_angstrom, template_to_fit, filts,
                    wv_corr, flux, time, errs, z, output_chi=True,
                    output_params=False))

    # Chooses the template that yields the lowest chi2
    gi = np.argmin(chi2)
    chi2 = chi2[gi]
    best_temp = templates[gi]

    return best_temp


def interpolate(lc, wv_corr, sn_type, use_mean, z, verbose, stepsize, kernel_width=None):
    '''
    Interpolate the LC using a 2D Gaussian Process (GP)

    Parameters
    ----------
    lc : numpy.array
        LC array
    wv_corr : float
        mean of wavelengths, needed to find wavelenght in angstroms
    sn_type : string
        type of supernova template being used for GP mean function
    use_mean : bool
        If True, use a template for GP mean function
        If False, use 0 as GP mean function
    z : float
        redshift
    stepsize : float
        Step size in days used for GP sampling
    kernel_width : list of float
        The width (:math:`r^2`) of the GP kernel in the (time, wavelength) direction.
        If not given, the kernel width will be optimized.

    Output
    ------
    dense_lc : numpy.array
        GP-interpolated LC and errors
    test_y : numpy.array
        A set of flux and wavelength values from the template, to be plotted
    test_times : numpy.array
        A set of time values to plot the template data against
    dense_times : numpy.array
        Highly sampled array of phases
    '''

    # Extract light curve parameters
    times, fluxes, wv_effs, errs, widths = lc

    stacked_data = np.vstack([times, wv_effs]).T
    ufilts = np.unique(wv_effs)
    ufilts_in_angstrom = ufilts*1000.0 + wv_corr
    nfilts = len(ufilts)
    dense_times = np.arange(int(np.floor(np.min(times))),
                           (int(np.ceil(np.max(times)))+1),stepsize)
    length_of_times = len(dense_times)
    x_pred = np.zeros((length_of_times*nfilts, 2))
    dense_fluxes = np.zeros((length_of_times, nfilts))
    dense_errs = np.zeros((length_of_times, nfilts))

    # test_y is only used if mean = True
    # but I still need it to exist either way
    test_y = []
    test_times = []
    if use_mean:
        template = generate_template(ufilts_in_angstrom, sn_type)
        if verbose:
            print('Fitting Template...')
        f_stretch, t_shift, t_stretch = fit_template(ufilts_in_angstrom,
                                                     template, wv_effs,
                                                     wv_corr, fluxes, times,
                                                     errs, z)

        # George needs the mean function to be in this format
        class snModel(Model):

            def get_value(self, param):
                t = (param[:, 0] * 1./t_stretch) + t_shift
                wv = param[:, 1]
                return np.asarray([template(*p)[0] for p in zip(t, wv)]) \
                    + f_stretch

        # Get Test data so that the template can be plotted
        mean = snModel()
        for i in ufilts_in_angstrom:
            test_wv = np.full((1, length_of_times), i)
            test_times = np.arange(int(np.floor(np.min(times))),
                                   int(np.ceil(np.max(times))+1), stepsize)
            test_x = np.vstack((test_times, test_wv)).T
            test_y.append(mean.get_value(test_x))
        test_y = np.asarray(test_y)

    # Set up gp
    kernel = np.var(fluxes) \
        * george.kernels.Matern32Kernel(kernel_width or [12., 0.1], ndim=2)
    if not use_mean:
        gp = george.GP(kernel, mean=0)
    else:
        gp = george.GP(kernel, mean=snModel())
    gp.compute(stacked_data, errs)

    if kernel_width is None:
        def neg_ln_like(p):
            gp.set_parameter_vector(p)
            return -gp.log_likelihood(fluxes)

        def grad_neg_ln_like(p):
            gp.set_parameter_vector(p)
            return -gp.grad_log_likelihood(fluxes)

        # Optimize gp parameters
        bnds = ((None, None), (None, None), (None, None))
        result = minimize(neg_ln_like,
                          gp.get_parameter_vector(),
                          jac=grad_neg_ln_like,
                          bounds = bnds)
        gp.set_parameter_vector(result.x)

    # Populate arrays with time and wavelength values to be fed into gp
    for jj, time in enumerate(np.arange(int(np.floor(np.min(times))),
                                        int(np.ceil(np.max(times)))+1, stepsize)):
        x_pred[jj*nfilts: jj*nfilts+nfilts, 0] = [time] * nfilts
        x_pred[jj*nfilts: jj*nfilts+nfilts, 1] = ufilts

    # Run gp to estimate interpolation
    pred, pred_var = gp.predict(fluxes, x_pred, return_var=True)

    # Populate dense_lc with newly gp-predicted values
    for jj in np.arange(nfilts):
        gind = np.where(np.abs(x_pred[:, 1]-ufilts[jj]) < epsilon)[0]
        dense_fluxes[:, int(jj)] = pred[gind]
        dense_errs[:, int(jj)] = np.sqrt(pred_var[gind])
    dense_lc = np.dstack((dense_fluxes, dense_errs))

    return dense_lc, test_y, test_times, dense_times


def fit_bb(dense_lc, wvs, use_mcmc, T_max):
    '''
    Fit a series of BBs to the GP LC
    Adapted from superbol, Nicholl, M. 2018, RNAAS)

    Parameters
    ----------
    dense_lc : numpy.array
        GP-interpolated LC
    wvs : numpy.array
        Reference wavelengths in Ang

    Output
    ------
    T_arr : numpy.array
        BB temperature array (K)
    R_arr : numpy.array
        BB radius array (cm)
    Terr_arr : numpy.array
        BB radius error array (K)
    Rerr_arr : numpy.array
        BB temperature error array (cm)
    '''

    T_arr = np.zeros(len(dense_lc))
    R_arr = np.zeros(len(dense_lc))
    Terr_arr = np.zeros(len(dense_lc))
    Rerr_arr = np.zeros(len(dense_lc))
    covar_arr = np.zeros(len(dense_lc))

    prior_fit = (9000, 1e15)

    for i, datapoint in enumerate(dense_lc):
        fnu = 10.**((-datapoint[:, 0]+48.6) / -2.5)
        ferr = datapoint[:, 1]
        fnu = fnu * 4. * np.pi * (3.086e19)**2
        fnu_err = np.abs(0.921034 * 10.**(0.4*datapoint[:, 0] - 19.44)) \
            * ferr * 4. * np.pi * (3.086e19)**2

        flam = fnu*c / (wvs*ang_to_cm)**2
        flam_err = fnu_err*c / (wvs*ang_to_cm)**2

        if use_mcmc:
            def log_likelihood(params, lam, f, f_err):
                T, R = params
                model = bbody(lam, T, R)
                return -np.sum((f-model)**2/(f_err**2))

            def log_prior(params):
                T, R = params
                if T > 0 and T < T_max and R > 0:
                    return 0.
                return -np.inf

            def log_probability(params, lam, f, f_err):
                lp = log_prior(params)
                if not np.isfinite(lp):
                    return -np.inf
                return lp + log_likelihood(params, lam, f, f_err)

            nwalkers = 16
            ndim = 2
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,
                                            args=[wvs, flam, flam_err])
            T0 = 9000 + 1000*np.random.rand(nwalkers)
            R0 = 1e15 + 1e14*np.random.rand(nwalkers)
            p0 = np.vstack([T0, R0])
            p0 = p0.T
            burn_in_state = sampler.run_mcmc(p0, 100)
            sampler.reset()
            sampler.run_mcmc(burn_in_state, 4000)
            flat_samples = sampler.get_chain(discard=100, thin=1, flat=True)
            covar_arr[i] = np.cov(flat_samples.T)[0,1]
            T_arr[i] = np.median(flat_samples[:, 0])
            R_arr[i] = np.median(flat_samples[:, 1])
            Terr_arr[i] = (np.percentile(flat_samples[:, 0], 84) -
                           np.percentile(flat_samples[:, 0], 16)) / 2.
            Rerr_arr[i] = (np.percentile(flat_samples[:, 1], 84) -
                           np.percentile(flat_samples[:, 1], 16)) / 2.

        else:

            try:
                BBparams, covar = curve_fit(bbody, wvs, flam, maxfev=10000,
                                            p0=prior_fit, sigma=flam_err,
                                            bounds=(0, [T_max, np.inf]),
                                            method = 'dogbox', 
                                            absolute_sigma=True)

                # Get temperature and radius, with errors, from fit
                T_arr[i] = BBparams[0]
                Terr_arr[i] = np.sqrt(np.diag(covar))[0]
                R_arr[i] = np.abs(BBparams[1])
                Rerr_arr[i] = np.sqrt(np.diag(covar))[1]
                covar_arr[i] = covar[0,1]
                prior_fit = BBparams
            except RuntimeWarning:
                T_arr[i] = np.nan
                R_arr[i] = np.nan
                Terr_arr[i] = np.nan
                Rerr_arr[i] = np.nan
                covar_arr[i] = np.nan

    return T_arr, R_arr, Terr_arr, Rerr_arr, covar_arr


def plot_gp(lc, dense_times, dense_lc, snname, flux_corr, my_filters, wvs, test_data,
            outdir, sn_type, test_times, mean, show_template):
    '''
    Plot the GP-interpolate LC and save

    Parameters
    ----------
    lc : numpy.array
        Original LC data
    dense_times : numpy.array
        Highly sampled array of phases
    dense_lc : numpy.array
        GP-interpolated LC
    snname : string
        SN Name
    flux_corr : float
        Flux correction factor for GP
    my_filters : list
        List of filters
    wvs : numpy.array
        List of central wavelengths, for colors
    outdir : string
        Output directory
    sn_type : string
        Type of sn template used for GP mean function
    test_times : numpy array
        Time values for sn template to be plotted against
    mean : bool
        Whether or not a non-zero mean function is being used in GP
    show_template : bool
        Whether or not the sn template is plotted

    Output
    ------
    '''

    # Extract light curve parameters
    times, fluxes, wv_effs, errs, widths = lc

    # Import a color map to make the plots look pretty
    cm = plt.get_cmap('rainbow')
    wv_colors = (wvs-np.min(wvs)) / (np.max(wvs)-np.min(wvs))

    # Plot interpolation, template, and error (shaded areas)
    for jj in np.arange(len(wv_colors)):
        plt.plot(dense_times, -dense_lc[:, jj, 0], color=cm(wv_colors[jj]),
                 label=my_filters[jj].split('/')[-1])
        if mean:
            if show_template:
                plt.plot(test_times, -(test_data[jj, :] + flux_corr), '--',
                         color=cm(wv_colors[jj]))  # Template curves
        plt.fill_between(dense_times,
                         -dense_lc[:, jj, 0] - dense_lc[:, jj, 1],
                         -dense_lc[:, jj, 0] + dense_lc[:, jj, 1],
                         color=cm(wv_colors[jj]), alpha=0.2)

    # Plot original data points and error bars
    for i, filt in enumerate(np.unique(wv_effs)):
        gind = np.where(wv_effs == filt)
        x = times[gind]
        x = x.flatten()
        y = -fluxes[gind] - flux_corr
        y = y.flatten()
        yerr = errs[gind]
        yerr = yerr.flatten()

        plt.plot(x, y, 'o', color=cm(wv_colors[i]))
        plt.errorbar(x, y, yerr=yerr, fmt='none', color=cm(wv_colors[i]))

    if mean:
        plt.title(snname + ' using sn' + sn_type + ' template')
    else:
        plt.title(snname + ' Light Curves')
    plt.legend()
    plt.xlabel('Time(days)')
    plt.ylabel('Absolute Magnitudes')
    plt.gca().invert_yaxis()
    plt.savefig(outdir + snname + '_' + str(sn_type) + '_gp.png')
    plt.clf()

    return 1


def plot_bb_ev(dense_times, Tarr, Rarr, Terr_arr, Rerr_arr, snname, outdir, sn_type):
    '''
    Plot the BB temperature and radius as a function of time

    Parameters
    ----------
    dense_times : numpy.array
        Highly sampled array of phases
    T_arr : numpy.array
        BB temperature array (K)
    R_arr : numpy.array
        BB radius array (cm)
    Terr_arr : numpy.array
        BB radius error array (K)
    Rerr_arr : numpy.array
        BB temperature error array (cm)
    snname : string
        SN Name
    outdir : string
        Output directory
    sn_type : string
        The type of sn template used for the gp

    Output
    ------
    '''

    fig, axarr = plt.subplots(2, 1, sharex=True)

    axarr[0].plot(dense_times, Tarr / 1.e3, color='k')
    axarr[0].fill_between(dense_times, Tarr/1.e3 - Terr_arr/1.e3,
                          Tarr/1.e3 + Terr_arr/1.e3, color='k', alpha=0.2)
    axarr[0].set_ylabel('Temp. (1000 K)')

    axarr[1].plot(dense_times, Rarr / 1e15, color='k')
    axarr[1].fill_between(dense_times, Rarr/1e15 - Rerr_arr/1e15,
                          Rarr/1e15 + Rerr_arr/1e15, color='k', alpha=0.2)
    axarr[1].set_ylabel(r'Radius ($10^{15}$ cm)')

    axarr[1].set_xlabel('Time (Days)')
    axarr[0].set_title(snname + ' Black Body Evolution')

    plt.savefig(outdir + snname + '_' + str(sn_type) + '_bb_ev.png')
    plt.clf()

    return 1


def plot_bb_bol(dense_times, bol_lum, bol_err, snname, outdir, sn_type):
    '''
    Plot the BB bolometric luminosity as a function of time

    Parameters
    ----------
    dense_times : numpy.array
        Highly sampled array of phases
    bol_lum : numpy.array
        BB bolometric luminosity (erg/s)
    bol_err : numpy.array
        BB bolometric luminosity error (erg/s)
    snname : string
        SN Name
    outdir : string
        Output directory
    sn_type : string
        The type of sn template used in the gp

    Output
    ------
    '''

    plt.plot(dense_times, bol_lum, color='k')
    plt.fill_between(dense_times, bol_lum-bol_err, bol_lum+bol_err,
                     color='k', alpha=0.2)

    plt.title(snname + ' Bolometric Luminosity')
    plt.xlabel('Time (Days)')
    plt.ylabel('Bolometric Luminosity')
    plt.yscale('log')
    plt.savefig(outdir + snname + '_' + str(sn_type) + '_bb_bol.png')
    plt.clf()

    return 1


def write_output(lc, dense_times, dense_lc, Tarr, Terr_arr, Rarr, Rerr_arr,
                 bol_lum, bol_err, ufilts,
                 snname, outdir, sn_type):
    '''
    Write out the interpolated LC and BB information

    Parameters
    ----------
    lc : numpy.array
        Initial light curve
    dense_times : numpy.array
        Highly sampled array of phases
    dense_lc : numpy.array
        GP-interpolated LC
    T_arr : numpy.array
        BB temperature array (K)
    Terr_arr : numpy.array
        BB radius error array (K)
    R_arr : numpy.array
        BB radius array (cm)
    Rerr_arr : numpy.array
        BB temperature error array (cm)
    bol_lum : numpy.array
        BB luminosity (erg/s)
    bol_err : numpy.array
        BB luminosity error (erg/s)
    ufilts : list
        List of filter names
    snname : string
        SN Name
    outdir : string
        Output directory
    sn_type : string
        Type of sn template used for the gp

    Output
    ------
    '''

    dense_lc = np.reshape(dense_lc, (len(dense_lc), -1))
    dense_lc = np.hstack((np.reshape(dense_times, (len(dense_times), 1)), dense_lc))
    tabledata = np.stack((Tarr, Terr_arr, Rarr, Rerr_arr, bol_lum, bol_err)).T
    tabledata = np.hstack((dense_lc, tabledata)).T
    table_header = ['Phase']
    for filt in ufilts:
        table_header.append(filt)
        table_header.append(filt + '_err')
    format_dict = {head: '%0.3f' for head in table_header}
    data_columns = ['Temp. (K)', 'Temp. Err. (K)',
                    'Radius (cm)', 'Radius Err. (cm)',
                    'Bol. Lum. (erg/s)', 'Bol. Err. (erg/s)']
    table_header.extend(data_columns)
    format_dict.update({head: '%0.3e' for head in data_columns})
    table = Table([*tabledata],
                  names=table_header,
                  meta={'name': 'first table'})
    for filt in ufilts:
        table[filt] *= -1.
    ascii.write(table, outdir + snname + '_' + str(sn_type) + '.txt',
                formats=format_dict, overwrite=True)

    return 1


def main(args=None):
    if args is None:
        args = sys.argv[1:]
    
    default_data = importlib_resources.files('extrabol.example') / 'SN2010bc.snana.dat'
    default_data = str(default_data)
    # Define all arguments
    parser = argparse.ArgumentParser(description='extrabol helpers')
    parser.add_argument('snfile', nargs='?',
                        default=default_data,
                        type=str, help='Give name of SN file')
    parser.add_argument('-m', '--mean', dest='mean', type=str, default='0',
                        help="Template function for gp.\
                                Choose '1a', '1bc', '2l', '2p', or \
                                '0' for no template")
    parser.add_argument('-t', '--show_template', dest='template',
                        action='store_true',
                        help="Shows template function on plots")
    parser.add_argument('-d', '--dist', dest='distance', type=float,
                        help='Object luminosity distance in Mpc')
    parser.add_argument('-z', '--redshift', dest='redshift', type=float,
                        help='Object redshift')
    parser.add_argument('--dm', dest='dm', type=float,
                        help='Object distance modulus')
    parser.add_argument("--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument("--plot", help="Make plots", dest='plot',
                        action="store_true", default=True)
    parser.add_argument("--outdir", help="Output directory", dest='outdir',
                        type=str, default='./products/')
    parser.add_argument("--ebv", help="Milky Way E(B-V)", dest='ebv',
                        type=float, default=0.)
    parser.add_argument("--hostebv", help="Host-galaxy E(B-V)", dest='hostebv',
                        type=float, default=0.)
    parser.add_argument('-s', '--start',
                        help='The time of the earliest \
                            data point to be accepted',
                        type=float, default=-50)
    parser.add_argument('-e', '--end',
                        help='The time of the latest \
                            data point to be accepted',
                        type=float, default=200)
    parser.add_argument('-sz', '--stepsize',
                        help='Step size in days of GP sampler',
                        type=float, default=0.1)
    parser.add_argument('-snr',
                        help='The minimum signal to \
                            noise ratio to be accepted',
                        type=float, default=4)
    parser.add_argument('-mc', '--use_mcmc', dest='mc',
                        help='Use a Markov Chain Monte Carlo \
                              to fit BBs instead of curve_fit. This will \
                              take longer, but gives better error estimates',
                        default=False, action="store_true")
    parser.add_argument('--T_max', dest='T_max',  help='Temperature prior \
                                                        for black body fits',
                        type=float, default=40000.)
    parser.add_argument('-k', '--kernel-width',
                        help='The width (:math:`r^2`) of the GP kernel in the (time, wavelength) direction. \
                              If not given, the kernel width will be optimized.',
                        type=float, nargs=2)

    args = parser.parse_args(args)

    # We need to know if an sn template is being used for gp
    sn_type = args.mean
    try:
        sn_type = int(sn_type)
        mean = False
        if sn_type != 0:
            print('Template request not valid. Assuming mean function of 0.')
    except ValueError:
        sn_type = sn_type
        mean = True

    # Read redshift and ebv from the file
    # *** EDIT HERE!!
    if '.snana' in args.snfile:
        metadata, data_df = read_snana_file(args.snfile)
        redshift = float(metadata['REDSHIFT'])
        ebv = float(metadata['MWEBV'])
    else:
        with open(args.snfile, 'r') as f:
            redshift = float(f.readline())
            ebv = float(f.readline())
    if args.redshift is not None:
        if args.redshift != redshift and args.verbose:
            print(f'Overriding redshift in file {redshift:f} with {args.redshift:f}')
        redshift = args.redshift
    if args.ebv is not None:
        if args.ebv != ebv and args.verbose:
            print(f'Overriding E(B-V) in file {ebv:f} with {args.ebv:f}')
        ebv = args.ebv

    # Solve for distance modulus if possible
    # if not, assume that data is already in absolute magnitudes
    if args.dm is not None:
        dm = args.dm
    elif args.distance is not None:
        dm = 5. * np.log10(args.distance * 1e5)
        if args.verbose:
            print(f'Calculating distance modulus {dm:f} from distance {args.distance:f} Mpc')
    elif redshift > 0.:
        dm = cosmo.distmod(redshift).value
        if args.verbose:
            print(f'Calculating distance modulus {dm:f} from redshift {redshift:f}')
    else:
        dm = 0.
        if args.verbose:
            print('Assuming absolute magnitudes.')

    # Make sure outdir name is formatted correctly
    if args.outdir[-1] != '/':
        args.outdir += '/'

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    snname = ('.').join(args.snfile.split('.')[: -1]).split('/')[-1]

    lc, wv_corr, flux_corr, my_filters = read_in_photometry(args.snfile,
                                                            dm,
                                                            redshift,
                                                            args.start,
                                                            args.end, args.snr,
                                                            ebv,
                                                            args.hostebv,
                                                            args.verbose)

    # Test which template fits the data best
    if sn_type == 'test':
        sn_type = test(lc, wv_corr, args.redshift)
    if args.verbose:
        print('Using ' + str(sn_type) + ' template.')

    dense_lc, test_data, test_times, dense_times = interpolate(lc, wv_corr, sn_type,
                                                  mean, args.redshift,
                                                  args.verbose, args.stepsize)
    wvs, wvind = np.unique(lc.T[:, 2], return_index=True)
    wvs = wvs*1000.0 + wv_corr
    my_filters = np.asarray(my_filters)
    ufilts = my_filters[wvind]

    # Converts to AB magnitudes
    dense_lc[:, :, 0] += flux_corr

    if args.verbose:
        print('Fitting Blackbodies, this may take a few minutes...')
    Tarr, Rarr, Terr_arr, Rerr_arr, covar_arr = fit_bb(dense_lc, wvs, args.mc,
                                                       args.T_max)

    # Calculate bolometric luminosity and error
    bol_lum = 4. * np.pi * Rarr**2 * sigsb * Tarr**4
    covar_err = 2. * (4. * np.pi * sigsb)**2 * (2 * Rarr * Tarr**4) * \
                (4 * Rarr**2 * Tarr**3) * covar_arr
    bol_err = 4. * np.pi * sigsb * np.sqrt(
                (2. * Rarr * Tarr**4 * Rerr_arr)**2
                + (4. * Tarr**3 * Rarr**2 * Terr_arr)**2
                )
    bol_err = np.sqrt(bol_err**2 + covar_err)

    if args.plot:
        if args.verbose:
            print('Making plots in ' + args.outdir)
        plot_gp(lc, dense_times, dense_lc, snname, flux_corr, ufilts, wvs, test_data,
                args.outdir, sn_type, test_times, mean, args.template)
        plot_bb_ev(dense_times, Tarr, Rarr, Terr_arr, Rerr_arr, snname,
                   args.outdir, sn_type)
        plot_bb_bol(dense_times, bol_lum, bol_err, snname, args.outdir, sn_type)

    if args.verbose:
        print('Writing output to ' + args.outdir)
    write_output(lc, dense_times, dense_lc, Tarr, Terr_arr, Rarr, Rerr_arr,
                 bol_lum, bol_err, ufilts, snname, args.outdir, sn_type)
    print('All done!')


if __name__ == "__main__":
    main()
