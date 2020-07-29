#!/usr/bin/env python

import numpy as np
from astroquery.svo_fps import SvoFps

def read_in_photometry(filename):
    '''
    Read in SN file
    '''
    photometry_data = np.loadtxt(filename,dtype=str)
    ufilts = np.unique(photometry_data[:,3])
    index = SvoFps.get_filter_index()
    filterIDs = np.asarray(index['filterID'].data,dtype=str)
    wavelengthEffs = np.asarray(index['WavelengthEff'].data,dtype=float)
    widthEffs = np.asarray(index['WidthEff'].data,dtype=float)
    zpts_all = np.asarray(index['ZeroPoint'].data,dtype=str)

    wv_effs = []
    width_effs = []
    for ufilt in ufilts:
        gind = np.where(filterIDs==ufilt)[0]
        wv_effs.append(wavelengthEffs[gind][0])
        width_effs.append(widthEffs[gind][0])    

    zpts = []
    for datapoint in photometry_data:
        if datapoint[-1] == 'AB':
            zpts.append(3631.00)
        else:
            gind = np.where(filterIDs==datapoint[3])
            zpts.append(float(zpts_all[gind[0]][0]))
    print(zpts)

    return 0

def interpolate(lc):
    '''
    Interpolate using 2D GP
    '''
    return 0


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
    read_in_photometry('./example/Gaia16apd.dat')


if __name__ == "__main__":
    main()