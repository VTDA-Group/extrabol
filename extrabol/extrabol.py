#!/usr/bin/env python

import numpy as np
from astroquery.svo_fps import SvoFps
import matplotlib.pyplot as plt
import george
from scipy.optimize import minimize, curve_fit
import argparse
from astropy.cosmology import Planck13 as cosmo
from astropy.cosmology import z_at_value
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
import pkg_resources

# Define a few important constants that will be used later
epsilon = 0.001
c = 2.99792458E10 # cm / s
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


def read_in_photometry(filename, dm, redshift, start, end, snr, mwebv, use_wc):
    '''
    Read in SN file

    Parameters
    ----------
    filename : string
        Input file name
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
        Extinction to be removed
    use_wc : bool
        If True, use redshift corrected wv for extinction correction

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

    photometry_data = np.loadtxt(filename, dtype=str, skiprows=2)

    # Extract key information into seperate arrays
    phases = np.asarray(photometry_data[:, 0], dtype=float)
    errs = np.asarray(photometry_data[:, 2], dtype=float)
    index = SvoFps.get_filter_index(timeout=3600)
    filterIDs = np.asarray(index['filterID'].data, dtype=str)
    wavelengthEffs = np.asarray(index['WavelengthEff'].data, dtype=float)
    widthEffs = np.asarray(index['WidthEff'].data, dtype=float)
    zpts_all = np.asarray(index['ZeroPoint'].data, dtype=str)

    # Extract filter names and effective wavelengths
    wv_effs = []
    width_effs = []
    my_filters = []
    for ufilt in photometry_data[:, 3]:
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
    for datapoint in photometry_data:
        mag = float(datapoint[1]) - dm
        if datapoint[-1] == 'AB':
            zpts.append(3631.00)
        else:
            gind = np.where(filterIDs == datapoint[3])
            zpts.append(float(zpts_all[gind[0]][0]))

        flux = 10.**(mag/-2.5) * zpts[-1] * (1.+redshift)

        # Convert Flux to log-flux space
        # This is easier on the Gaussian Process
        # 'fluxes' is also equivilant to the negative absolute magnitude
        flux = 2.5 * (np.log10(flux)-np.log10(3631.00))
        fluxes.append(flux)

    # Remove extinction
    if use_wc:
        ext = extinction.fm07(wv_effs / (1.+redshift), mwebv)
    else:
        ext = extinction.fm07(wv_effs, mwebv)
    for i in np.arange(len(fluxes)):
        fluxes[i] = fluxes[i] + ext[i]

    # The GP prefers when values are relatively close to zero
    # so we adjust wavelength and flux accordingly
    # This will be undone after interpolation
    wv_corr = np.mean(wv_effs / (1.+redshift))
    flux_corr = np.min(fluxes) - 1.0
    wv_effs = wv_effs - wv_corr
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
    phases = np.asarray(phases) - phases[peak_i]

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

    my_template_file = pkg_resources.resource_filename('extrabol.template_bank', 'smoothed_sn' + sn_type + '.npz')
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

    #set peak flux to t=0
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


def interpolate(lc, wv_corr, sn_type, use_mean, z):
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

    Output
    ------
    dense_lc : numpy.array
        GP-interpolated LC and errors
    test_y : numpy.array
        A set of flux and wavelength values from the template, to be plotted
    test_times : numpy.array
        A set of time values to plot the template data against
    '''

    lc = lc.T

    times = lc[:, 0]
    fluxes = lc[:, 1]
    wv_effs = lc[:, 2]
    errs = lc[:, 3]
    stacked_data = np.vstack([times, wv_effs]).T
    ufilts = np.unique(lc[:, 2])
    ufilts_in_angstrom = ufilts*1000.0 + wv_corr
    nfilts = len(ufilts)
    x_pred = np.zeros((int((np.ceil(np.max(times))-np.floor(np.min(times)))+1)*nfilts, 2))
    dense_fluxes = np.zeros((int((np.ceil(np.max(times))-np.floor(np.min(times)))+1), nfilts))
    dense_errs = np.zeros((int((np.ceil(np.max(times))-np.floor(np.min(times)))+1), nfilts))

    # test_y is only used if mean = True
    # but I still need it to exist either way
    test_y = []
    test_times = []
    if use_mean:
        template = generate_template(ufilts_in_angstrom, sn_type)
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
            test_wv = np.full((1, round(np.max(times))
                               - round(np.min(times))), i)
            test_times = np.arange(round(np.min(times)), round(np.max(times)))
            test_x = np.vstack((test_times, test_wv)).T
            test_y.append(mean.get_value(test_x))
        test_y = np.asarray(test_y)

    # Set up gp
    kernel = np.var(lc[:, 1]) \
        * george.kernels.ExpSquaredKernel([50, 0.5], ndim=2)
    if not use_mean:
        gp = george.GP(kernel, mean=0)
    else:
        gp = george.GP(kernel, mean=snModel())
    gp.compute(stacked_data, lc[:, -2])

    def neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.log_likelihood(lc[:, 1])

    def grad_neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.grad_log_likelihood(lc[:, 1])

    # Optomize gp parameters
    result = minimize(neg_ln_like,
                      gp.get_parameter_vector(),
                      jac=grad_neg_ln_like)
    gp.set_parameter_vector(result.x)

    # Populate arrays with time and wavelength values to be fed into gp
    for jj, time in enumerate(np.arange(int(np.floor(np.min(times))), int(np.ceil(np.max(times)))+1)):
        x_pred[jj*nfilts: jj*nfilts+nfilts, 0] = [time] * nfilts
        x_pred[jj*nfilts: jj*nfilts+nfilts, 1] = ufilts

    # Run gp to estimate interpolation
    pred, pred_var = gp.predict(lc[:, 1], x_pred, return_var=True)

    # Populate dense_lc with newly gp-predicted values
    for jj in np.arange(nfilts):
        gind = np.where(np.abs(x_pred[:, 1]-ufilts[jj]) < epsilon)[0]
        dense_fluxes[:, int(jj)] = pred[gind]
        dense_errs[:, int(jj)] = np.sqrt(pred_var[gind])
    dense_lc = np.dstack((dense_fluxes, dense_errs))

    return dense_lc, test_y, test_times


def fit_bb(dense_lc, wvs, use_mcmc):
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
                if T > 0 and T < 20000. and R > 0:
                    return 0.
                return -np.inf
            def log_probability(params, lam, f, f_err):
                lp = log_prior(params)
                if not np.isfinite(lp):
                    return -np.inf
                return lp + log_likelihood(params, lam, f, f_err)

            nwalkers = 16
            ndim = 2
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=[wvs,flam,flam_err])
            T0 = 9000 + 1000*np.random.rand(nwalkers)
            R0 = 1e15 + 1e14*np.random.rand(nwalkers)
            p0 = np.vstack([T0,R0])
            p0 = p0.T

            burn_in_state = sampler.run_mcmc(p0, 100)
            sampler.reset()
            sampler.run_mcmc(burn_in_state, 4000)
            flat_samples = sampler.get_chain(discard=100, thin=1, flat=True)

            T_arr[i] = np.median(flat_samples[:,0])
            R_arr[i] = np.median(flat_samples[:,1])
            Terr_arr[i] = (np.percentile(flat_samples[:,0], 84)-np.percentile(flat_samples[:,0], 16)) / 2.
            Rerr_arr[i] = (np.percentile(flat_samples[:,1], 84)-np.percentile(flat_samples[:,1], 16)) / 2.

        else:
            try:
                BBparams, covar = curve_fit(bbody, wvs, flam, maxfev=8000,
                                            p0=(9000, 1e15), sigma=flam_err,
                                            bounds=(0, [20000, np.inf]))
                # Get temperature and radius, with errors, from fit
                T_arr[i] = BBparams[0]
                Terr_arr[i] = np.sqrt(np.diag(covar))[0]
                R_arr[i] = np.abs(BBparams[1])
                Rerr_arr[i] = np.sqrt(np.diag(covar))[1]
            except RuntimeWarning:
                T_arr[i] = np.nan
                R_arr[i] = np.nan
                Terr_arr[i] = np.nan
                Rerr_arr[i] = np.nan


    return T_arr, R_arr, Terr_arr, Rerr_arr


def plot_gp(lc, dense_lc, snname, flux_corr, my_filters, wvs, test_data,
            outdir, sn_type, test_times, mean, show_template):
    '''
    Plot the GP-interpolate LC and save

    Parameters
    ----------
    lc : numpy.array
        Original LC data
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
    plot_times = np.arange(int(np.floor(np.min(lc[:,0]))), (int(np.ceil(np.max(lc[:,0])))+1))

    # Import a color map to make the plots look pretty
    cm = plt.get_cmap('rainbow')
    wv_colors = (wvs-np.min(wvs)) / (np.max(wvs)-np.min(wvs))

    # Plot interpolation, template, and error (shaded areas)
    for jj in np.arange(len(wv_colors)):
        plt.plot(plot_times, -dense_lc[:, jj, 0], color=cm(wv_colors[jj]),
                 label=my_filters[jj].split('/')[-1])
        if mean:
            if show_template:
                plt.plot(test_times, -(test_data[jj, :] + flux_corr), '--',
                         color=cm(wv_colors[jj]))  # Template curves
        plt.fill_between(plot_times,
                         -dense_lc[:, jj, 0] - dense_lc[:, jj, 1],
                         -dense_lc[:, jj, 0] + dense_lc[:, jj, 1],
                         color=cm(wv_colors[jj]), alpha=0.2)

    # Plot original data points and error bars
    for i, filt in enumerate(np.unique(lc[:, 2])):
        gind = np.where(lc[:, 2] == filt)
        x = lc[gind, 0]
        x = x.flatten()
        y = -lc[gind, 1] - flux_corr
        y = y.flatten()
        yerr = lc[gind, 3]
        yerr = yerr.flatten()

        plt.plot(x, y, 'o', color=cm(wv_colors[i]))
        plt.errorbar(x, y, yerr=yerr, fmt='none', color=cm(wv_colors[i]))

    if mean:
        plt.title(snname + ' using sn' + sn_type + ' template')
    else:
        plt.title(snname)
    plt.legend()
    plt.xlabel('Time(days)')
    plt.ylabel('Absolute Magnitudes')
    plt.gca().invert_yaxis()
    plt.savefig(outdir + snname + '_' + str(sn_type) + '_gp.png')
    plt.clf()

    return 1


def plot_bb_ev(lc, Tarr, Rarr, Terr_arr, Rerr_arr, snname, outdir, sn_type):
    '''
    Plot the BB temperature and radius as a function of time

    Parameters
    ----------
    lc : numpy.array
        Original LC data
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
    interp_times = np.arange(int(np.floor(np.min(lc[:,0]))), int(np.ceil(np.max(lc[:,0])))+1)
    fig, axarr = plt.subplots(2, 1, sharex=True)

    axarr[0].plot(interp_times, Tarr / 1.e3, color='k')
    axarr[0].fill_between(interp_times, Tarr/1.e3 - Terr_arr/1.e3, Tarr/1.e3 + Terr_arr/1.e3, color='k', alpha=0.2)
    axarr[0].set_ylabel('Temp. (1000 K)')

    axarr[1].plot(interp_times, Rarr / 1e15, color='k')
    axarr[1].fill_between(interp_times, Rarr/1e15 - Rerr_arr/1e15, Rarr/1e15 + Rerr_arr/1e15, color='k', alpha=0.2)
    axarr[1].set_ylabel(r'Radius ($10^{15}$ cm)')

    axarr[1].set_xlabel('Time (Days)')
    axarr[0].set_title(snname)

    plt.savefig(outdir + snname + '_' + str(sn_type) + '_bb_ev.png')
    plt.clf()

    return 1


def plot_bb_bol(lc, bol_lum, bol_err, snname, outdir, sn_type):
    '''
    Plot the BB bolometric luminosity as a function of time

    Parameters
    ----------
    lc : numpy.array
        Original LC data
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
    plot_times = np.arange(int(np.floor(np.min(lc[:,0]))), int(np.ceil(np.max(lc[:,0])))+1)

    plt.plot(plot_times, bol_lum, color='k')
    plt.fill_between(plot_times, bol_lum-bol_err, bol_lum+bol_err, color='k', alpha=0.2)

    plt.title(snname)
    plt.xlabel('Time (Days)')
    plt.ylabel('Bolometric Luminosity')
    plt.yscale('log')
    plt.savefig(outdir + snname + '_' + str(sn_type) + '_bb_bol.png')
    plt.clf()

    return 1


def write_output(lc, dense_lc, Tarr, Terr_arr, Rarr, Rerr_arr,
                 bol_lum, bol_err, my_filters,
                 snname, outdir, sn_type):
    '''
    Write out the interpolated LC and BB information

    Parameters
    ----------
    lc : numpy.array
        Initial light curve
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
    my_filters : list
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

    times = np.arange(int(np.floor(np.min(lc[:,0]))), int(np.ceil(np.max(lc[:,0])))+1)
    dense_lc = np.reshape(dense_lc, (len(dense_lc), -1))
    dense_lc = np.hstack((np.reshape(-times, (len(times), 1)), dense_lc))
    tabledata = np.stack((Tarr / 1e3, Terr_arr / 1e3, Rarr / 1e15,
                          Rerr_arr / 1e15, np.log10(bol_lum),
                          np.log10(bol_err))).T
    tabledata = np.hstack((-dense_lc, tabledata)).T

    ufilts = np.unique(my_filters)
    table_header = []
    table_header.append('Time (MJD)')
    for filt in ufilts:
        table_header.append(filt)
        table_header.append(filt + '_err')
    table_header.extend(['Temp./1e3 (K)', 'Temp. Err.',
                         'Radius/1e15 (cm)', 'Radius Err.',
                         'Log10(Bol. Lum)', 'Log10(Bol. Err)'])
    table = Table([*tabledata],
                  names=table_header,
                  meta={'name': 'first table'})

    format_dict = {head: '%0.3f' for head in table_header}
    ascii.write(table, outdir + snname + '_' + str(sn_type) + '.txt', formats=format_dict,
                overwrite=True)

    return 1


def main():

    default_data = pkg_resources.resource_filename('extrabol.example', 'PSc000174_extrabol.dat')
    # Define all arguments
    parser = argparse.ArgumentParser(description='extrabol helpers')
    parser.add_argument('snfile', nargs='?',
                        default=default_data,
                        type=str, help='Give name of SN file')
    parser.add_argument('-m', '--mean', dest='mean', type=str, default='0',
                        help="Template function for gp.\
                                Choose \'1a\',\'1bc\', \'2l\', \'2p\', or \
                                \'0\' for no template")
    parser.add_argument('-t', '--show_template', dest='template',
                        action='store_true',
                        help="Shows template function on plots")
    parser.add_argument('-d', '--dist', dest='distance', type=float,
                        help='Object luminosity distance', default=1e-5)
    parser.add_argument('-z', '--redshift', dest='redshift', type=float,
                        help='Object redshift', default=-1.)
                        # Redshift can't =-1
                        # this is simply a flag to be replaced later
    parser.add_argument('--dm', dest='dm', type=float, default=0,
                        help='Object distance modulus')
    parser.add_argument("--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument("--plot", help="Make plots", dest='plot',
                        action="store_true", default=True)
    parser.add_argument("--outdir", help="Output directory", dest='outdir',
                        type=str, default='./products/')
    parser.add_argument("--ebv", help="MWebv", dest='ebv',
                        type=float, default=-1.)
                        # Ebv won't =-1
                        # this is another flag to be replaced later
    parser.add_argument("--hostebv", help="Host B-V", dest='hostebv',
                        type=float, default=0.0)
    parser.add_argument('-s', '--start',
                        help='The time of the earliest \
                            data point to be accepted',
                        type=float, default=-50)
    parser.add_argument('-e', '--end',
                        help='The time of the latest \
                            data point to be accepted',
                        type=float, default=200)
    parser.add_argument('-snr',
                        help='The minimum signal to \
                            noise ratio to be accepted',
                        type=float, default=4)
    parser.add_argument('-wc', '--wvcorr',
                        help='Use the redshift-corrected \
                            wavelenghts for extinction calculations',
                        action="store_true")
    parser.add_argument('-mc', '--use_mcmc', dest='mc',
                        help='Use a Markov Chain Monte Carlo \
                              to fit BBs instead of curve_fit. This will \
                              take longer, but gives better error estimates',
                        default=False, action="store_true")

    args = parser.parse_args()

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

    # If redshift or ebv aren't specified by the user,
    # we read them in from the file here
    if args.redshift == -1 or args.ebv == -1:
        # Read in redshift and ebv and replace values if not specified
        f = open(args.snfile, 'r')
        if args.redshift == -1:
            args.redshift = float(f.readline())
            if args.ebv == -1:
                args.ebv = float(f.readline())
        if args.ebv == -1:
            args.ebv = float(f.readline())
            args.ebv = float(f.readline())
        f.close

    # Solve for redshift, distance, and/or dm if possible
    # if not, assume that data is already in absolute magnitudes
    if args.redshift != 0 or args.distance != 1e-5 or args.dm != 0:
        if args.redshift != 0:
            args.distance = cosmo.luminosity_distance(args.redshift).value
            args.dm = cosmo.distmod(args.redshift).value
        elif args.distance != 1e-5:
            args.redshift = z_at_value(cosmo.luminosity_distance, distance
                                       * u.Mpc)
            dm = cosmo.distmod(args.redshift).value
        else:
            args.redshift = z_at_value(cosmo.distmod, dm * u.mag)
            distance = cosmo.luminosity_distance(args.redshift).value
    elif args.verbose:
        print('Assuming absolute magnitudes.')

    # Make sure outdir name is formatted correctly
    if args.outdir[-1] != '/':
        args.outdir += '/'

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    snname = ('.').join(args.snfile.split('.')[: -1]).split('/')[-1]

    lc, wv_corr, flux_corr, my_filters = read_in_photometry(args.snfile,
                                                            args.dm,
                                                            args.redshift,
                                                            args.start,
                                                            args.end, args.snr,
                                                            args.ebv,
                                                            args.wvcorr)

    # Test which template fits the data best
    if sn_type == 'test':
        sn_type = test(lc, wv_corr, args.redshift)
    if args.verbose:
        print('Using ' + str(sn_type) + ' template.')

    dense_lc, test_data, test_times = interpolate(lc, wv_corr, sn_type,
                                                  mean, args.redshift)
    lc = lc.T

    wvs, wvind = np.unique(lc[:, 2], return_index=True)
    wvs = wvs*1000.0 + wv_corr
    my_filters = np.asarray(my_filters)
    ufilts = my_filters[wvind]

    # Converts to AB magnitudes
    dense_lc[:, :, 0] += flux_corr

    Tarr, Rarr, Terr_arr, Rerr_arr = fit_bb(dense_lc, wvs, args.mc)

    # Calculate bolometric luminosity and error
    bol_lum = 4. * np.pi * Rarr**2 * sigsb * Tarr**4
    bol_err = 4. * np.pi * sigsb * np.sqrt(
                (2. * Rarr * Tarr**4 * Rerr_arr)**2
                + (4. * Tarr**3 * Rarr**2 * Terr_arr)**2
                )

    if args.plot:
        if args.verbose:
            print('Making plots in ' + args.outdir)
        plot_gp(lc, dense_lc, snname, flux_corr, ufilts, wvs, test_data,
                args.outdir, sn_type, test_times, mean, args.template)
        plot_bb_ev(lc, Tarr, Rarr, Terr_arr, Rerr_arr, snname, args.outdir, sn_type)
        plot_bb_bol(lc, bol_lum, bol_err, snname, args.outdir, sn_type)

    if args.verbose:
        print('Writing output to ' + args.outdir)
    write_output(lc, dense_lc, Tarr, Terr_arr, Rarr, Rerr_arr,
                 bol_lum, bol_err, my_filters, snname, args.outdir, sn_type)
    print('job completed')


if __name__ == "__main__":
    main()