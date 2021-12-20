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

epsilon = 0.0001
c = 2.99792458E10
sigsb = 5.6704e-5 #erg / cm^2 / s / K^4
h = 6.62607E-27
ang_to_cm = 1e-8
k_B = 1.38064852E-16 # cm^2 * g / s^2 / K

def bbody(lam,T,R):
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
    exponential = (h * c) / (lam_cm * k_B * T)
    blam = ((2. * np.pi * h * c ** 2) / (lam_cm ** 5)) / (np.exp(exponential) - 1.)
    area = 4. * np.pi * R**2
    lum = blam * area

    return lum

def read_in_photometry(filename, dm, z):
    '''
    Read in SN file

    Parameters
    ----------
    filename : string
        Input file name
    dm : float
        Distance modulus
    z : float
        redshift

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
    photometry_data = np.loadtxt(filename,dtype=str)

    phases = np.asarray(photometry_data[:,0],dtype=float)
    errs = np.asarray(photometry_data[:,2],dtype=float)
    index = SvoFps.get_filter_index(timeout=3600)
    filterIDs = np.asarray(index['filterID'].data,dtype=str)
    wavelengthEffs = np.asarray(index['WavelengthEff'].data,dtype=float)
    widthEffs = np.asarray(index['WidthEff'].data,dtype=float)
    zpts_all = np.asarray(index['ZeroPoint'].data,dtype=str)

    wv_effs = []
    width_effs = []
    my_filters = []
    for ufilt in photometry_data[:,3]:
        gind = np.where(filterIDs==ufilt)[0]
        if len(gind) == 0:
            sys.exit('Cannot find '+str(ufilt)+' in SVO.')
        wv_effs.append(wavelengthEffs[gind][0])
        width_effs.append(widthEffs[gind][0])   
        my_filters.append(ufilt) 

    zpts = []
    fluxes = []
    for datapoint in photometry_data:
        mag = float(datapoint[1]) - dm
        if datapoint[-1] == 'AB':
            zpts.append(3631.00)
        else:
            gind = np.where(filterIDs==datapoint[3])
            zpts.append(float(zpts_all[gind[0]][0]))
    
        flux = 10.** (mag / -2.5) * zpts[-1] * (1.+z)
        #I'm calling ths flux..but I'm going to do it in log flux space
        flux = 2.5 * (np.log10(flux) - np.log10(3631.00))
        fluxes.append(flux)

    wv_effs = np.asarray(wv_effs)

    wv_corr = np.mean(wv_effs/(1.+z))
    flux_corr = np.min(fluxes) - 1.0
    wv_effs = wv_effs - wv_corr
    fluxes = np.asarray(fluxes) - flux_corr
    phases = np.asarray(phases) - np.min(phases)
    lc = np.vstack((phases,fluxes,wv_effs/1000.,errs,width_effs))

    return lc,wv_corr,flux_corr, my_filters

def interpolate(lc):
    '''
    Interpolate the LC using a 2D Gaussian Process (GP)

    Parameters
    ----------
    lc : numpy.array
        LC array

    Output
    ------
    dense_lc : numpy.array
        GP-interpolated LC and errors
    '''
    lc = lc.T

    template = np.load('./smoothed_sn1a.npz')
    temp_times = template['time']
    temp_wavelength = template['wavelength']
    temp_f_lambda = template['f_lambda']

    #template_function = interp.interp2d(temp_times, temp_wavelength, temp_f_lambda)

    times = lc[:,0]
    filters = lc[:,2]
    stacked_data = np.vstack([times, filters]).T
    ufilts = np.unique(lc[:,2])
    nfilts = len(ufilts)
    x_pred = np.zeros((len(lc)*nfilts, 2))
    dense_fluxes = np.zeros((len(lc), nfilts))
    dense_errs = np.zeros((len(lc), nfilts))
    kernel = np.var(lc[:,1]) * george.kernels.ExpSquaredKernel([50, 0.5], ndim=2)
    gp = george.GP(kernel, mean = 0)
    gp.compute(stacked_data, lc[:,-2])

    def neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.log_likelihood(lc[:,1])

    def grad_neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.grad_log_likelihood(lc[:,1])

    result = minimize(neg_ln_like,
                      gp.get_parameter_vector(),
                      jac=grad_neg_ln_like)
    gp.set_parameter_vector(result.x)

    for jj, time in enumerate(lc[:,0]):
        x_pred[jj*nfilts:jj*nfilts+nfilts, 0] = [time]*nfilts
        x_pred[jj*nfilts:jj*nfilts+nfilts, 1] = ufilts

    pred, pred_var = gp.predict(lc[:,1], x_pred, return_var=True)

    for jj in np.arange(nfilts):
        gind = np.where(np.abs(x_pred[:, 1]-ufilts[jj])<epsilon)[0]
        dense_fluxes[:, int(jj)] = pred[gind]
        dense_errs[:, int(jj)] = np.sqrt(pred_var[gind])
    dense_lc = np.dstack((dense_fluxes, dense_errs))
    
    return dense_lc


def fit_bb(dense_lc,wvs):
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

    for i,datapoint in enumerate(dense_lc):

        fnu = 10.**((-datapoint[:,0] + 48.6)/-2.5)
        ferr = datapoint[:,1]
        fnu = fnu * 4. * np.pi * (3.086e19)**2 #LAZZYYYY assumption of 10 kpc
        fnu_err = np.abs(0.921034 * 10.**(0.4 * datapoint[:,0] - 19.44)) * \
                ferr * 4. * np.pi * (3.086e19)**2

        flam = fnu * c / (wvs * ang_to_cm)**2 
        flam_err = fnu_err * c / (wvs * ang_to_cm)**2

        try:
            BBparams, covar = curve_fit(bbody,wvs,flam,
                                p0=(9000,1e15),sigma=flam_err)
            # Get temperature and radius, with errors, from fit
            T1 = BBparams[0]
            T1_err = np.sqrt(np.diag(covar))[0]
            R1 = np.abs(BBparams[1])
            R1_err = np.sqrt(np.diag(covar))[1]
        except:
            T1 = np.nan
            R1 = np.nan
            T1_err = np.nan
            R1_err = np.nan

        T_arr[i] = T1
        R_arr[i] = R1
        Terr_arr[i] = T1_err
        Rerr_arr[i] = R1_err

    return T_arr,R_arr,Terr_arr,Rerr_arr

def plot_gp(lc, dense_lc, snname, flux_corr, my_filters, wvs, outdir):
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

    Output
    ------
    '''

    cm = plt.get_cmap('rainbow') 
    wv_colors = (wvs - np.min(wvs)) / (np.max(wvs) - np.min(wvs))

    gind = np.argsort(lc[:,0])
    for jj in np.arange(len(wv_colors)):
        plt.plot(lc[gind,0],-dense_lc[gind,jj,0],color=cm(wv_colors[jj]),
                label=my_filters[jj].split('/')[-1])
        plt.fill_between(lc[gind,0],-dense_lc[gind,jj,0]-dense_lc[gind,jj,1],
                    -dense_lc[gind,jj,0]+dense_lc[gind,jj,1],
                    color=cm(wv_colors[jj]),alpha=0.2)

    for i,filt in enumerate(np.unique(lc[:,2])):
        gind = np.where(lc[:,2]==filt)
        plt.plot(lc[gind,0],-(lc[gind,1] + flux_corr),'o',
                color=cm(wv_colors[i]))
    plt.title(snname)
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Absolute Magnitudes')
    # Uhg, magnitudes are the worst.
    plt.gca().invert_yaxis()
    plt.savefig(outdir+snname+'_gp.png')
    plt.clf()
    return 1

def plot_bb_ev(lc, Tarr, Rarr, Terr_arr, Rerr_arr, snname, outdir):
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

    Output
    ------
    '''
    fig,axarr = plt.subplots(2,1,sharex=True)
    axarr[0].plot(lc[:,0],Tarr/1.e3,'ko')
    axarr[0].errorbar(lc[:,0],Tarr/1.e3,yerr=Terr_arr/1.e3,fmt='none',color='k')
    axarr[0].set_ylabel('Temp. (1000 K)')

    axarr[1].plot(lc[:,0],Rarr/1e15,'ko')
    axarr[1].errorbar(lc[:,0],Rarr/1e15,yerr=Rerr_arr/1e15,fmt='none',color='k')
    axarr[1].set_ylabel(r'Radius ($10^{15}$ cm)')

    axarr[1].set_xlabel('Time (Days)')
    axarr[0].set_title(snname)

    plt.savefig(outdir+snname+'_bb_ev.png')
    plt.clf()
    return 1

def plot_bb_bol(lc, bol_lum, bol_err, snname, outdir):
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

    Output
    ------
    '''
    plt.plot(lc[:,0],bol_lum,'ko')
    plt.errorbar(lc[:,0],bol_lum,yerr=bol_err,fmt='none',color='k')

    plt.title(snname)
    plt.xlabel('Time (Days)')
    plt.ylabel('Bolometric Luminosity')
    plt.yscale('log')
    plt.savefig(outdir+snname+'_bb_bol.png')
    plt.clf()
    return 1

def write_output(lc, dense_lc,Tarr,Terr_arr,Rarr,Rerr_arr,
                 bol_lum,bol_err,my_filters,
                 snname,outdir):
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

    Output
    ------
    '''
    times = lc[:,0]
    dense_lc = np.reshape(dense_lc,(len(dense_lc),-1))
    dense_lc = np.hstack((np.reshape(-times,(len(times),1)),dense_lc))
    tabledata = np.stack((Tarr/1e3,Terr_arr/1e3,Rarr/1e15,
                        Rerr_arr/1e15, np.log10(bol_lum), np.log10(bol_err))).T
    tabledata = np.hstack((-dense_lc,tabledata)).T

    ufilts = np.unique(my_filters)
    table_header = []
    table_header.append('Time (MJD)')
    for filt in ufilts:
        table_header.append(filt)
        table_header.append(filt+'_err')
    table_header.extend(['Temp./1e3 (K)','Temp. Err.',
        'Radius/1e15 (cm)','Radius Err.', 'Log10(Bol. Lum)', 'Log10(Bol. Err)'])
    table = Table([*tabledata],
        names = table_header,
                meta={'name': 'first table'})

    format_dict = {head:'%0.3f' for head in table_header}
    ascii.write(table, outdir+snname+'.txt', formats=format_dict, overwrite=True)
    return 1

def main(snfile, distance=1e-5, redshift = 0, dm = 0, verbose=True, plot=True, outdir = 'outdir', 
                ebv = 0.0, hostebv = 0.0):

    if (redshift != 0) | (distance != 1e-5) | (dm != 0):
        if redshift !=0 :
            distance = cosmo.luminosity_distance(redshift).value
            dm = cosmo.distmod(redshift).value
        elif distance != 1e-5:
            redshift = z_at_value(cosmo.luminosity_distance,distance * u.Mpc)
            dm = cosmo.distmod(redshift).value
        else:
            redshift = z_at_value(cosmo.distmod,dm * u.mag)
            distance = cosmo.luminosity_distance(redshift).value
    elif verbose:
        print('Assuming absolute magnitudes.')

    if outdir[-1] != '/':
        outdir+='/'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    snname = ('.').join(snfile.split('.')[:-1]).split('/')[-1]

    lc,wv_corr,flux_corr, my_filters = read_in_photometry(snfile, dm, redshift)
    dense_lc = interpolate(lc)
    lc = lc.T

    wvs,wvind = np.unique(lc[:,2],return_index=True)
    wvs = wvs * 1000.0 + wv_corr
    my_filters = np.asarray(my_filters)
    ufilts = my_filters[wvind]


    dense_lc[:,:,0] += flux_corr # This is now in AB mags

    Tarr,Rarr,Terr_arr,Rerr_arr = fit_bb(dense_lc,wvs)

    bol_lum = 4. * np.pi * Rarr **2 * sigsb * Tarr**4
    bol_err = 4. * np.pi * sigsb * np.sqrt(\
                (2. * Rarr * Tarr**4 * Rerr_arr)**2 + \
                (4. * Tarr**3 * Rarr**2 * Terr_arr)**2)

    if plot:
        if verbose:
            print('Making plots in '+outdir)
        plot_gp(lc,dense_lc,snname,flux_corr,ufilts,wvs,outdir)
        plot_bb_ev(lc,Tarr,Rarr,Terr_arr,Rerr_arr,snname,outdir)
        plot_bb_bol(lc, bol_lum, bol_err, snname, outdir)

    if verbose:
        print('Writing output to '+outdir)
    write_output(lc,dense_lc,Tarr,Terr_arr,Rarr,Rerr_arr,bol_lum,bol_err,my_filters,
                snname,outdir)


if __name__ == "__main__":
    main('./extrabol/example/Gaia16apd.dat',  dm = 38.38)
