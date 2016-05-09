"""
This fucntion write out the parameters that were used for a cloudy_slug run.
"""

import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")

def write_cluster_cloudyparams(data, model_name, fmt):
    """
    Write out photometry for nebular emission computed by cloudy on a
    slug spectrum for a series of clusters

    Parameters
       data : namedtuple
          Cluster cloudy parameter data; a namedtuple containing the
          fields id, trial, time, cloudy_hden, cloudy_r0, cloudy_rS, 
          cloudy_QH0, cloudy_covFac, cloudy_U, and cloudy_Omega
       model_name : string
          Base file name to give the model to be written. Can include a
          directory specification if desired.
       fmt : string
          Format for the output file. Allowed values are 'ascii', 'bin'
          or 'binary, and 'fits'.

    Returns
       Nothing
    """

    # Make sure fmt is valid
    if fmt != 'ascii' and fmt != 'bin' and fmt != 'binary' and \
       fmt != 'fits':
        raise ValueError("fmt must be ascii, bin, binary, or fits")

    # Make sure we're not trying to do fits if we don't have astropy
    if fmt == 'fits' and fits is None:
        raise ValueError("Couldn't import astropy, so fits format "+
                         "is unavailable.")

    if fmt == 'ascii':

        # ASCII mode
        fp = open(model_name+'_cluster_cloudyparams.txt', 'w')

        # Write header lines
        fp.write(("{:<14s}"*9).
                 format('UniqueID', 'Time', 'hden', 'r0', 'rS',
                        'QH0', 'CovFac', 'U', 'Omega') + "\n")
        fp.write(("{:<14s}"*9).
                 format('', '(yr)', '(cm^-3)', '(cm)', '(s^-1)',
                        '', '', '', '') + "\n")
        fp.write(("{:<14s}"*9).
                 format('-----------', '-----------', '-----------',
                        '-----------', '-----------', '-----------',
                        '-----------', '-----------', '-----------')
                 + "\n")

        # Write data
        for i in range(data.trial.size):
            # If this is a new trial, write a separator
            if i != 0:
                if data.trial[i] != data.trial[i-1]:
                    fp.write("-"*(9*14-3)+"\n")
            fp.write(("{:11d}   {:11.5e}   {:11.5e}   {:11.5e}   " +
                      "{:11.5e}   {:11.5e}   {:11.5e}   " +
                      "{:11.5e}   {:11.5e}\n").format(
                          data.id[i], data.time[i],
                          data.cloudy_hden[i], data.cloudy_r0[i],
                          data.cloudy_rS[i],
                          data.cloudy_QH0[i], data.cloudy_covFac[i],
                          data.cloudy_U[i], data.cloudy_Omega[i]))

        # Close
        fp.close()

        
    elif fmt == 'bin' or fmt == 'binary':

        # Binary mode
        fp = open(model_name+'_cluster_cloudyparams.bin', 'wb')

        # Break data into blocks of clusters with the same time
        # and trial number
        ptr = 0
        while ptr < data.trial.size:

            # Find the next cluster that differs from this one in
            # either time or trial number
            diff = np.where(
                np.logical_or(data.trial[ptr+1:] != data.trial[ptr],
                              data.time[ptr+1:] != data.time[ptr]))[0]
            if diff.size == 0:
                block_end = data.trial.size
            else:
                block_end = ptr + diff[0] + 1

            # Write out time and number of clusters
            ncluster = block_end - ptr
            fp.write(np.uint(data.trial[ptr]))
            fp.write(data.time[ptr])
            fp.write(ncluster)

            # Loop over clusters and write them
            for k in range(ptr, block_end):
                fp.write(data.cloudy_hden[k])
                fp.write(data.cloudy_r0[k])
                fp.write(data.cloudy_rS[k])
                fp.write(data.cloudy_QH0[k])
                fp.write(data.cloudy_covFac[k])
                fp.write(data.cloudy_U[k])
                fp.write(data.cloudy_Omega[k])
        
            # Move pointer
            ptr = block_end
                
        # Close file
        fp.close()
   
    elif fmt == 'fits' or fmt == 'fits2':

        # Convert data to FITS columns
        cols = []
        cols.append(fits.Column(name="Trial", format="1K",
                                unit="", array=data.trial))
        cols.append(fits.Column(name="UniqueID", format="1K",
                                unit="", array=data.id))
        cols.append(fits.Column(name="Time", format="1D",
                                unit="yr", array=data.time))
        cols.append(fits.Column(name="hden", format="1D",
                                unit="cm^-3", array=data.cloudy_hden))
        cols.append(fits.Column(name="r0", format="1D",
                                unit="cm", array=data.cloudy_r0))
        cols.append(fits.Column(name="rS", format="1D",
                                unit="cm", array=data.cloudy_rS))
        cols.append(fits.Column(name="QH0",
                                format="1D",
                                unit="s^-1", array=data.cloudy_QH0))
        cols.append(fits.Column(name="covFac", format="1D",
                                unit="", array=data.cloudy_covFac))
        cols.append(fits.Column(name="U", format="1D",
                                unit="", array=data.cloudy_U))
        cols.append(fits.Column(name="Omega", format="1D",
                                unit="", array=data.cloudy_Omega))
        fitscols = fits.ColDefs(cols)

        # Create the binary table HDU
        tbhdu = fits.BinTableHDU.from_columns(fitscols)

        # Create dummy primary HDU
        prihdu = fits.PrimaryHDU()

        # Create HDU list and write to file
        hdulist = fits.HDUList([prihdu, tbhdu])
        hdulist.writeto(model_name+'_cluster_cloudyparams.fits',
                        clobber=True)

