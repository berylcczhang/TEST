"""
This fucntion writes out the parameters that were used for a cloudy_slug run.
"""

import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")

def write_integrated_cloudyparams(data, model_name, fmt):
    """
    Write out photometry for nebular emission computed by cloudy on a
    slug spectrum for a series of clusters

    Parameters
       data : namedtuple
          Cluster cloudy parameter data; a namedtuple containing the
          fields time, cloudy_hden, cloudy_r0, cloudy_rS, 
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
        fp = open(model_name+'_integrated_cloudyparams.txt', 'w')

        # Write header lines
        fp.write(("{:<14s}"*8).
                 format('Time', 'Hden', 'R0', 'RS',
                        'QH0', 'CovFac', 'U', 'Omega') + "\n")
        fp.write(("{:<14s}"*8).
                 format('(yr)', '(cm^-3)', '(cm)', '(s^-1)',
                        '', '', '', '') + "\n")
        fp.write(("{:<14s}"*8).
                 format('-----------', '-----------', '-----------',
                        '-----------', '-----------', '-----------',
                        '-----------', '-----------')
                 + "\n")
        
        # Write data
        ntime = data.cloudy_hden.shape[-2]
        ntrial = data.cloudy_hden.shape[-1]
        if len(data.time) > ntime:
            random_time = True
        else:
            random_time = False
        for i in range(ntrial):
            if i != 0:
                fp.write("-"*(8*14-3)+"\n")
            for j in range(ntime):
                if random_time:
                    t_out = data.time[i]
                else:
                    t_out = data.time[j]
                fp.write(("{:11.5e}   {:11.5e}   {:11.5e}   " +
                          "{:11.5e}   {:11.5e}   {:11.5e}   " +
                          "{:11.5e}   {:11.5e}\n")
                         .format(t_out, data.cloudy_hden[j,i],
                                 data.cloudy_r0[j,i],
                                 data.cloudy_rS[j,i],
                                 data.cloudy_QH0[j,i],
                                 data.cloudy_covFac[j,i],
                                 data.cloudy_U[j,i],
                                 data.cloudy_Omega[j,i]))

        # Close
        fp.close()
        
    elif fmt == 'bin' or fmt == 'binary':

        # Binary mode
        fp = open(model_name+'_integrated_cloudyparams.bin', 'wb')

        # Write data
        ntime = data.cloudy_hden.shape[-2]
        ntrial = data.cloudy_hden.shape[-1]
        if len(data.time) > ntime:
            random_time = True
        else:
            random_time = False
        for i in range(ntrial):
            for j in range(ntime):
                if random_time:
                    t_out = data.time[i]
                else:
                    t_out = data.time[j]
                fp.write(np.uint(i))
                fp.write(t_out)
                fp.write(data.cloudy_hden[j,i])
                fp.write(data.cloudy_r0[j,i])
                fp.write(data.cloudy_rS[j,i])
                fp.write(data.cloudy_QH0[j,i])
                fp.write(data.cloudy_covFac[j,i])
                fp.write(data.cloudy_U[j,i])
                fp.write(data.cloudy_Omega[j,i])

        # Close
        fp.close()

    elif fmt == 'fits':

        # FITS mode

        # Figure out number of trials, and tile arrays
        ntrial = data.cloudy_hden.shape[-1]
        ntimes = data.cloudy_hden.shape[-2]
        trial = np.transpose(np.tile(
            np.arange(ntrial, dtype='int64'), (ntimes,1))).\
            flatten()
        if len(data.time) > ntimes:
            times = data.time
        else:
            times = np.tile(data.time, ntrial)

        # Convert data to FITS columns
        cols = []
        cols.append(fits.Column(name="Trial", format="1K",
                                unit="", array=trial))
        cols.append(fits.Column(name="Time", format="1D",
                                unit="yr", array=times))
        cols.append(fits.Column(name="HDen", format="1D",
                                unit="cm^-3",
                                array=np.transpose(data.cloudy_hden).
                                flatten()))
        cols.append(fits.Column(name="R0", format="1D",
                                unit="cm",
                                array=np.transpose(data.cloudy_r0).
                                flatten()))
        cols.append(fits.Column(name="RS", format="1D",
                                unit="cm",
                                array=np.transpose(data.cloudy_rS).
                                flatten()))
        cols.append(fits.Column(name="QH0", format="1D",
                                unit="s^-1",
                                array=np.transpose(data.cloudy_QH0).
                                flatten()))
        cols.append(fits.Column(name="covFac", format="1D",
                                array=np.transpose(data.cloudy_covFac).
                                flatten()))
        cols.append(fits.Column(name="U", format="1D",
                                array=np.transpose(data.cloudy_U).
                                flatten()))
        cols.append(fits.Column(name="Omega", format="1D",
                                array=np.transpose(data.cloudy_Omega).
                                flatten()))
        fitscols = fits.ColDefs(cols)

        # Create the binary table HDU
        tbhdu = fits.BinTableHDU.from_columns(fitscols)

        # Create dummy primary HDU
        prihdu = fits.PrimaryHDU()

        # Create HDU list and write to file
        hdulist = fits.HDUList([prihdu, tbhdu])
        hdulist.writeto(model_name+'_integrated_cloudyparams.fits',
                        clobber=True)
