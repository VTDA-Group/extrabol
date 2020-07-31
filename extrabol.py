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
from astropy.table import QTable
from astropy.io import ascii

epsilon = 0.0001
c = 2.99792458E10
sigsb = 5.6704e-5 #erg / cm^2 / s / K^4

# STOLEN FROM MATT
def bbody(lam,T,R):
    '''
    Calculate the corresponding blackbody radiance for a set
    of wavelengths given a temperature and radiance.
    Parameters
    ---------------
    lam: Reference wavelengths in Angstroms
    T:   Temperature in Kelvin
    R:   Radius in cm
    Output
    ---------------
    Spectral radiance in units of erg/s/Angstrom
    (calculation and constants checked by Sebastian Gomez)
    '''

    #R = R * 1e15 #scale for code
    #T = T * 1e4 #scale for code

    # Planck Constant in cm^2 * g / s
    h = 6.62607E-27


    # Convert wavelength to cm
    lam_cm = lam * 1E-8

    # Boltzmann Constant in cm^2 * g / s^2 / K
    k_B = 1.38064852E-16

    # Calculate Radiance B_lam, in units of (erg / s) / cm ^ 2 / cm
    exponential = (h * c) / (lam_cm * k_B * T)
    B_lam = ((2. * np.pi * h * c ** 2) / (lam_cm ** 5)) / (np.exp(exponential) - 1.)

    # Multiply by the surface area
    A = 4. * np.pi * R**2

    # Output radiance in units of (erg / s) / cm
    Radiance = B_lam * A

    return Radiance

def read_in_photometry(filename, dm, z):
    '''
    Read in SN file
    '''
    photometry_data = np.loadtxt(filename,dtype=str)

    phases = np.asarray(photometry_data[:,0],dtype=float)
    errs = np.asarray(photometry_data[:,2],dtype=float)
    index = SvoFps.get_filter_index()
    filterIDs = np.asarray(index['filterID'].data,dtype=str)
    wavelengthEffs = np.asarray(index['WavelengthEff'].data,dtype=float)
    widthEffs = np.asarray(index['WidthEff'].data,dtype=float)
    zpts_all = np.asarray(index['ZeroPoint'].data,dtype=str)

    wv_effs = []
    width_effs = []
    my_filters = []
    for ufilt in photometry_data[:,3]:
        gind = np.where(filterIDs==ufilt)[0]
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


    wv_corr = np.mean(wv_effs/(1.+z))
    flux_corr = np.min(fluxes) - 1.0
    wv_effs = np.asarray(wv_effs) - wv_corr
    fluxes = np.asarray(fluxes) - flux_corr
    phases = np.asarray(phases)
    lc = np.vstack((phases,fluxes,wv_effs/1000.,errs,width_effs))

    return lc,wv_corr,flux_corr, my_filters

def interpolate(lc):
    '''
    Interpolate using 2D GP
    '''
    lc = lc.T

    times = lc[:,0]
    filters = lc[:,2]
    stacked_data = np.vstack([times, filters]).T#get everything but width eff..hacky
    ufilts = np.unique(lc[:,2])
    nfilts = len(ufilts)
    x_pred = np.zeros((len(lc)*nfilts, 2))
    dense_fluxes = np.zeros((len(lc), nfilts))
    dense_errs = np.zeros((len(lc), nfilts))
    kernel = np.var(lc[:,1]) * george.kernels.ExpSquaredKernel([10, 0.5], ndim=2)
    gp = george.GP(kernel)
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
    Fit each interpolate point to a BB, taking into account GP error
    '''

    #STOLEN FROM MATT!!!
    # Fit blackbody to SED (the one that is not padded with zeros)

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

        flam = fnu * c / (wvs * 1e-8)**2 
        flam_err = fnu_err * c / (wvs * 1e-8)**2

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

def plot_gp(lc, dense_lc, snname, flux_corr, outdir):
    '''
    Make nice outputs
    '''

    colors = ['blue','red','green','yellow','purple','cyan']
    gind = np.argsort(lc[:,0])
    for jj in np.arange(6):
        plt.plot(lc[gind,0],dense_lc[gind,jj,0],color=colors[jj])
        plt.fill_between(lc[gind,0],dense_lc[gind,jj,0]-dense_lc[gind,jj,1],
                    dense_lc[gind,jj,0]+dense_lc[gind,jj,1],
                    color=colors[jj],alpha=0.2)

    for i,filt in enumerate(np.unique(lc[:,2])):
        gind = np.where(lc[:,2]==filt)
        plt.plot(lc[gind,0],lc[gind,1] + flux_corr,'o',color=colors[i])
    plt.title(snname)
    plt.xlabel('Time')
    plt.ylabel('Negative Absolute Magnitudes')
    plt.savefig(outdir+snname+'_gp.png')
    plt.clf()
    return 1

def plot_bb_ev(lc, Tarr, Rarr, Terr_arr, Rerr_arr, snname, outdir):
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
    plt.plot(lc[:,0],np.log10(bol_lum),'ko')
    #plt.errorbar(lc[:,0],bol_lum,yerr=bol_err,fmt='none',color='k')

    plt.title(snname)
    plt.xlabel('Time (Days)')
    plt.ylabel('Bolometric Luminosity')
    #plt.yscale('log')
    plt.savefig(outdir+snname+'_bb_bol.png')
    plt.clf()
    return 1

def main():
    parser = argparse.ArgumentParser(description='extrabol helpers')
    parser.add_argument('snfile', type=str, help='Give name of SN file')
    parser.add_argument('-d','--dist', dest='distance', type=float,
                    help='Object luminosity distance', default=1e-5)
    parser.add_argument('-z','--redshift', dest='redshift', type=float,
                    help='Object redshift', default = 0)
    parser.add_argument('--dm', dest='dm', type=float, default=0,
                    help='Object distance modulus')
    parser.add_argument("--verbose", help="increase output verbosity",
                    action="store_true")
    parser.add_argument("--plot", help="Make plots",dest='plot',
                    type=bool, default=True)
    parser.add_argument("--outdir", help="Output directory",dest='outdir',
                    type=str, default='./products/')
    parser.add_argument("--ebv", help="MWebv",dest='ebv',
                    type=float, default=0.0)
    parser.add_argument("--hostebv", help="Host B-V",dest='hostebv',
                    type=float, default=0.0)
    
    args = parser.parse_args()

    if (args.redshift != 0) | (args.distance != 1e-5) | (args.dm != 0):
        if args.redshift !=0 :
            args.distance = cosmo.luminosity_distance(args.redshift).value
            args.dm = cosmo.distmod(args.redshift).value
        elif args.distance != 1e-5:
            args.redshift = z_at_value(cosmo.luminosity_distance,args.distance * u.Mpc)
            args.dm = cosmo.distmod(args.redshift).value
        else:
            args.redshift = z_at_value(cosmo.distmod,args.dm * u.mag)
            args.distance = cosmo.luminosity_distance(args.redshift).value
    elif args.verbose:
        print('Assuming absolute magnitudes.')

    if args.outdir[-1] != '/':
        args.outdir+='/'

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    snname = ('.').join(args.snfile.split('.')[:-1]).split('/')[-1]

    lc,wv_corr,flux_corr, my_filters = read_in_photometry(args.snfile, args.dm, args.redshift)
    dense_lc = interpolate(lc)
    lc = lc.T
    wvs = np.unique(lc[:,2]) * 1000.0 + wv_corr


    dense_lc[:,:,0] += flux_corr # This is now in AB mags

    Tarr,Rarr,Terr_arr,Rerr_arr = fit_bb(dense_lc,wvs)

    bol_lum = 4. * np.pi * Rarr **2 * sigsb * Tarr**4
    bol_err = 4. * np.pi * sigsb * np.sqrt(\
                (2. * Rarr * Tarr**4 * Rerr_arr)**2 + \
                (4. * Tarr**3 * Rarr**2 * Terr_arr)**2)

    if args.plot:
        if args.verbose:
            print('Making plots in '+args.outdir)
        plot_gp(lc,dense_lc,snname,flux_corr,args.outdir)
        plot_bb_ev(lc,Tarr,Rarr,Terr_arr,Rerr_arr,snname,args.outdir)
        plot_bb_bol(lc, bol_lum, bol_err, snname, args.outdir)

    #Final printout
    dense_lc = np.reshape(dense_lc,(len(dense_lc),-1))
    tabledata = np.stack((Tarr/1e3,Terr_arr/1e3,Rarr/1e15,Rerr_arr/1e15)).T
    tabledata = np.hstack((-dense_lc,tabledata)).T

    ufilts = np.unique(my_filters)
    table_header = []
    for filt in ufilts:
        table_header.append(filt)
        table_header.append(filt+'_err')
    table_header.extend(['Temp./1e3 (K)','Temp. Err.','Radius/1e15 (cm)','Radius Err.'])
    table = QTable([*tabledata],
        names = {(*table_header)},
                meta={'name': 'first table'})
    import sys
    format_dict = {head:'%0.3f' for head in table_header}
    ascii.write(table, args.outdir+snname+'.txt', formats=format_dict)
if __name__ == "__main__":
    main()