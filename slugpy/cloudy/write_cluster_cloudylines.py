"""
This function writes out cluster line luminosities computed by cloudy
from a slug run.
"""

from collections import namedtuple
import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")

def write_cluster_cloudylines(data, model_name, fmt):
    """
    Write out data computed by cloudy on a slug spectrum

    Parameters
       data : namedtuple
          Cloudy spectral data for clusters to be written; a namedtuple
          containing the fields time, cloudy_linelist, cloudy_linewl, 
          cloudy_linelum
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
        fp = open(model_name+'_cluster_cloudylines.txt', 'w')

        # Write header lines
        fp.write(("{:<14s}"*5).
                 format('UniqueID', 'Time', 'LineLabel', 'Wavelength', 
                        'Luminosity') + "\n")
        fp.write(("{:<14s}"*5).
                 format('', '(yr)', '', '(Angstrom)', '(erg/s)')
                 + "\n")
        fp.write(("{:<14s}"*5).
                 format('-----------', '-----------', '-----------',
                        '-----------', '-----------')
                 + "\n")

        # Write data
        for i in range(data.cloudy_linelum.shape[0]):
            # If this is a new trial, write a separator
            if i != 0:
                if data.trial[i] != data.trial[i-1]:
                    fp.write("-"*(5*14-3)+"\n")
            for j in range(data.cloudy_linelum.shape[1]):
                fp.write(("{:11d}   {:11.5e}   {:>11s}   {:11.5e}   " +
                          "{:11.5e}\n")
                         .format(data.id[i], data.time[i], 
                                 data.cloudy_linelist[j],
                                 data.cloudy_linewl[j],
                                 data.cloudy_linelum[i,j]))

        # Close
        fp.close()

    elif fmt == 'bin' or fmt == 'binary':

        # Binary mode
        fp = open(model_name+'_cluster_cloudylines.bin', 'wb')

        # Write out number of lines and line labels as ASCII, one per
        # line
        nlabel = len(data.cloudy_linelist)
        fp.write(str(nlabel)+"\n")
        for i in range(nlabel):
            fp.write(data.cloudy_linelist[i]+"\n")

        # Write line wavelengths
        fp.write(data.cloudy_linewl)

        # Construct lists of times and trials
        trials = np.unique(data.trial)
        times = np.unique(data.time)

        # Break data into blocks of clusters with the same time
        # and trial number
        ptr = 0
        for i in range(len(trials)):
            for j in range(len(times)):

                # Find block of clusters with this time and trial
                block_end = ptr + np.argmin(
                    np.logical_and(data.trial[ptr:] == trials[i],
                                   data.time[ptr:] == times[j]))

                # Special case: if block_end is the last entry in
                # the data, check if the last entry is the same as
                # the previous one. If so, move block_end one
                # space, to off the edge of the data.
                if block_end == len(data.trial)-1 and \
                   data.trial[-1] == trials[i] and \
                   data.time[-1] == times[j]:
                    block_end = block_end+1

                # Write out time and number of clusters
                ncluster = block_end - ptr
                fp.write(times[j])
                fp.write(ncluster)

                # Loop over clusters and write them
                for k in range(ptr, block_end):
                    fp.write(data.id[k])
                    fp.write(data.cloudy_linelum[k,:])

                # Move pointer
                ptr = block_end

        # Close file
        fp.close()

    elif fmt == 'fits':

        # FITS mode

        # Create a first HDU containing the line wavelengths and labels
        nl = len(data.cloudy_linewl)
        fmtstring = "A4"
        wlcols = [fits.Column(name="Line_Label",
                              format=fmtstring,
                              array=data.cloudy_linelist)]
        fmtstring = "D"
        wlcols.append(fits.Column(name="Wavelength",
                                  format=fmtstring,
                                  unit="Angstrom", 
                                  array=data.cloudy_linewl))
        wlfits = fits.ColDefs(wlcols)
        wlhdu = fits.BinTableHDU.from_columns(wlcols)

        # Convert line data to FITS columns, and make an HDU from
        # them
        speccols = []
        speccols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=data.trial))
        speccols.append(fits.Column(name="UniqueID", format="1K",
                                    unit="", array=data.id))
        speccols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=data.time))
        speccols.append(fits.Column(name="Line_Luminosity",
                                    format=str(nl)+"D",
                                    unit="erg/s",
                                    array=data.cloudy_linelum))
        specfits = fits.ColDefs(speccols)
        spechdu = fits.BinTableHDU.from_columns(specfits)

        # Create dummy primary HDU
        prihdu = fits.PrimaryHDU()

        # Create HDU list and write to file
        hdulist = fits.HDUList([prihdu, wlhdu, spechdu])
        hdulist.writeto(model_name+'_cluster_cloudylines.fits', 
                        clobber=True)
