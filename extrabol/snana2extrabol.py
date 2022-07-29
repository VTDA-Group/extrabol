
import numpy as np
import os
import argparse

def convert():
    parser = argparse.ArgumentParser(description='in and out directories')
    parser.add_argument('snana_in',
                        help='Specify the location of snana files to be converted')
    parser.add_argument('extrabol_out', 
                        help='Specify the location that extrabol files will be written to',
                        default='./extrabol_in/')
    args = parser.parse_args()

    directory = os.fsencode(args.snana_in)
    count = 0
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        snana = np.genfromtxt(args.snana_in+filename,
                                dtype=str, skip_header=17, skip_footer=1, usecols=(1,2,4,5))
        '''
        MJD     filt(griz)      FLUXCAL     FLUXCALERR
        '''
        mjd = np.asarray(snana[:,0], dtype=float)
        griz = np.asarray(snana[:,1], dtype=str)
        fluxcal = np.asarray(snana[:,2], dtype=float)
        fluxcalerr = np.asarray(snana[:,3], dtype=float)

        #Get redshift and mwebv

        f = open(args.snana_in+filename, 'r')
        for x in f:
            if f.readline(17) == 'REDSHIFT_FINAL:  ':
                redshift = f.readline(6)
        f.close
        f = open(args.snana_in+filename, 'r')
        for x in f:
            if f.readline(8) == 'MWEBV:  ':
                mwebv = f.readline(6)
        f.close

        #eliminate bad data

        gis=[]
        for i in np.arange(len(fluxcal)):
            if fluxcal[i] > 0:                          #flux needs to be positive
                gis.append(i)
        mjd = mjd[gis]
        griz = griz[gis]
        fluxcal = fluxcal[gis]
        fluxcalerr = fluxcalerr[gis]

        #time to build up the input for extrabol

        mag = -2.5 * np.log10(fluxcal) + 27.5
        magerr = fluxcalerr/fluxcal
        filt = []
        for i in np.arange(len(griz)):
            filt_t = 'PAN-STARRS/PS1.' + griz[i]
            filt.append(filt_t)
        type = np.full(len(mjd), 'AB')

        extrabol = np.vstack((mjd, mag, magerr, filt, type))
        extrabol = extrabol.T
        sn_name = filename.replace('PS1_PS1MD_', '')
        sn_name = sn_name.replace('.snana.dat', '')
        header = str(redshift) +'\n'+ str(mwebv)

        if not os.path.exists(args.extrabol_out):
            os.makedirs(args.extrabol_out)

        np.savetxt(args.extrabol_out + sn_name + '_extrabol.dat', extrabol, fmt = '%s', header = header, comments = '')
        count = count +1
        if count%500 == 0:
            print(str(count) + ' out of 5243')

        return 0

if __name__ == "__convert__":
    convert()