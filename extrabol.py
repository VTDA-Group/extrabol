#!/usr/bin/env python

import numpy as np
from astroquery.svo_fps import SvoFps
import matplotlib.pyplot as plt
import george
from scipy.optimize import minimize

epsilon = 0.0001

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
        mag = float(datapoint[1])
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


def fit_bb(dense_lc):
    '''
    Fit each interpolate point to a BB, taking into account GP error
    '''
    return 0

def make_outputs(bb_param):
    '''
    Make nice outputs
    '''




def main():
    lc,wv_corr,flux_corr = read_in_photometry('./example/Gaia16apd.dat')
    dense_lc = interpolate(lc)
    lc = lc.T

    dense_lc[:,:,0] += flux_corr # This is now in AB mags

    #This code doesn't work, just test
    colors = ['blue','red','green','yellow']
    gind = np.argsort(lc[:,0])
    for jj in np.arange(4):
        plt.plot(lc[gind,0],dense_lc[gind,jj,0],color=colors[jj])
        plt.fill_between(lc[gind,0],dense_lc[gind,jj,0]-dense_lc[gind,jj,1],
                    dense_lc[gind,jj,0]+dense_lc[gind,jj,1],
                    color=colors[jj],alpha=0.2)

    plt.plot(lc[gind,0],lc[gind,1] + flux_corr,'ko')
    plt.title('Gaia16apd')
    plt.xlabel('MJD')
    plt.ylabel('Negative Mag')
    plt.show()

if __name__ == "__main__":
    main()