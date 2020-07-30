#!/usr/bin/env python

import numpy as np
from astroquery.svo_fps import SvoFps
import matplotlib.pyplot as plt
import george
from scipy.optimize import minimize, curve_fit

epsilon = 0.0001
DM = 38.38 # HARD CODING CAUSE LAZY
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

def read_in_photometry(filename):
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
    for ufilt in photometry_data[:,3]:
        gind = np.where(filterIDs==ufilt)[0]
        wv_effs.append(wavelengthEffs[gind][0])
        width_effs.append(widthEffs[gind][0])    

    zpts = []
    fluxes = []
    for datapoint in photometry_data:
        mag = float(datapoint[1]) - DM
        if datapoint[-1] == 'AB':
            zpts.append(3631.00)
        else:
            gind = np.where(filterIDs==datapoint[3])
            zpts.append(float(zpts_all[gind[0]][0]))
    
        flux = 10.** (mag / -2.5) * zpts[-1]
        #I'm calling ths flux..but I'm going to do it in log flux space
        flux = 2.5 * (np.log10(flux) - np.log10(3631.00))
        fluxes.append(flux)


    wv_corr = np.mean(wv_effs)
    flux_corr = np.min(fluxes) - 1.0
    wv_effs = np.asarray(wv_effs) - wv_corr
    fluxes = np.asarray(fluxes) - flux_corr
    phases = np.asarray(phases)
    lc = np.vstack((phases,fluxes,wv_effs/1000.,errs,width_effs))

    return lc,wv_corr,flux_corr

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
    for i,datapoint in enumerate(dense_lc):

        fnu = 10.**((-datapoint[:,0] + 48.6)/-2.5)
        fnuplus = 10.**((-datapoint[:,0] + datapoint[:,1] + 48.6)/-2.5) 
        fnu = fnu * 4. * np.pi * (3.086e19)**2 #LAZZYYYY assumption of 10 kpc
        fnuplus = fnuplus * 4. * np.pi * (3.086e19)**2

        flam = fnu * c / (wvs * 1e-8)**2 
        flamplus = fnuplus * c / (wvs * 1e-8)**2
        ferr = np.abs(flam-flamplus)

        flam = flam 
        ferr = ferr / 1e87

        try:
            BBparams, covar = curve_fit(bbody,wvs,flam,
                                p0=(9000,1e15),sigma=ferr)
            # Get temperature and radius, with errors, from fit
            T1 = BBparams[0]
            T1_err = np.sqrt(np.diag(covar))[0]
            R1 = np.abs(BBparams[1])
            R1_err = np.sqrt(np.diag(covar))[1]
        except:
            T1 = np.nan
            R1 = np.nan

        
        
        plt.plot(wvs,flam,'k')
        plt.plot(wvs,bbody(wvs,T1,R1),'r')
        print(T1,R1)
        


        T_arr[i] = T1
        R_arr[i] = R1
    plt.show()
    return T_arr,R_arr

def make_outputs(bb_param):
    '''
    Make nice outputs
    '''


def main():
    lc,wv_corr,flux_corr = read_in_photometry('./example/Gaia16apd.dat')
    dense_lc = interpolate(lc)
    lc = lc.T

    wvs = np.unique(lc[:,2]) * 1000.0 + wv_corr


    dense_lc[:,:,0] += flux_corr # This is now in AB mags

    #This code doesn't work, just test
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
    plt.title('Gaia16apd')
    plt.xlabel('MJD')
    plt.ylabel('Negative Mag')
    plt.show()

    Tarr,Rarr = fit_bb(dense_lc,wvs)
    plt.plot(lc[:,0],Tarr,'o')
    plt.show()
    plt.plot(lc[:,0],Rarr,'o')
    plt.show()

    bol_lum = 4. * np.pi * Rarr **2 * sigsb * Tarr**4
    plt.plot(lc[:,0],bol_lum,'o')
    plt.title('Gaia16apd')
    plt.xlabel('MJD')
    plt.ylabel('Bolometric Luminosity')
    plt.show()


if __name__ == "__main__":
    main()