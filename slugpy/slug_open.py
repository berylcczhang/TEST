"""
Helper function to open slug output files.
"""

import os
import os.path as osp
import errno
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")

def slug_open(filename, output_dir=None, fmt=None):
    """
    Function to open a SLUG2 output file.

    Parameters
       filename : string
          Name of the file to open, without any extension. The following
          extensions are tried, in order: .txt, .bin, .fits
       output_dir : string
          The directory where the SLUG2 output is located; if set to None,
          the current directory is searched, followed by the
          SLUG_DIR/output directory if the SLUG_DIR environment variable
          is set
       fmt : 'txt' | 'ascii' | 'bin' | 'binary' | 'fits' | 'fits2'
          Format for the file to be read. If one of these is set, the
          function will only attempt to open ASCII-('txt' or 'ascii'), 
          binary ('bin' or 'binary'), or FITS ('fits' or 'fits2')
          formatted output, ending in .txt., .bin, or .fits,
          respectively. If set to None, the code will try to open
          ASCII files first, then if it fails try binary files, and if
          it fails again try FITS files.

    Returns
       fp : file or astropy.io.fits.hdu.hdulist.HDUList
          A file object pointing the file that has been opened
       fname : string
          Name of the file that was opened

    Raises
       IOError, if a file of the specified name cannot be found
    """

    # Did we get a specific directory in which to look? If not, try
    # current directory
    if output_dir is None:
        outdir = "."
    else:
        outdir = output_dir

    # Make sure fmt is valid
    if fmt != 'ascii' and fmt != 'txt' and fmt != 'bin' and \
       fmt != 'binary' and fmt != 'fits' and fmt != 'fits2' and \
              fmt is not None:
        raise ValueError("unknown format {}".fmt)

    # Make sure we're not trying to do fits if we don't have astropy
    if (fmt == 'fits' or fmt == 'fits2') and fits is None:
        raise ValueError("Couldn't import astropy, so fits format "+
                         "is unavailable.")

    # See if we have a text file
    if fmt is None or fmt=='ascii' or fmt=='txt':
        fname = osp.join(outdir, filename+'.txt')
        try:
            fp = open(fname, 'r')
            fmt = 'ascii'
        except IOError:
            fp = None
    else:
        fp = None

    # If that failed, look for a binary file
    if fp is None:
        if fmt is None or fmt=='bin' or fmt=='binary':
            fname = osp.join(outdir, filename+'.bin')
            try:
                fp = open(fname, 'rb')
                fmt = 'bin'
            except IOError:
                pass

    # If that failed, look for a fits file
    if fp is None:
        if fmt is None or fmt=='fits' or fmt=='fits2':
            fname = osp.join(outdir, filename+'.fits')
            try:
                fp = fits.open(fname)
                fmt = 'fits'
            except IOError:
                pass


    # If that failed, and we didn't get an explicit directory
    # specification, try looking in SLUG_DIR/output
    if (fp is None) and (output_dir is None) and \
       ('SLUG_DIR' in os.environ):
        outdir = osp.join(os.environ['SLUG_DIR'], 'output')
        if fmt is None or fmt == 'ascii' or fmt == 'txt':
            fname = osp.join(outdir,
                             filename+'.txt')
            try:
                fp = open(fname, 'r')
                fmt = 'ascii'
            except IOError:
                pass
        if (fmt is None or fmt == 'bin' or fmt == 'binary') \
           and fp is None:
            fname = osp.join(outdir, 
                             filename+'.bin')
            try:
                fp = open(fname, 'rb')
                fmt = 'bin'
            except IOError:
                pass
        if (fmt is None or fmt == 'fits' or fmt == 'fits2') \
           and fp is None and fits is not None:
            fname = osp.join(outdir, 
                             filename+'.fits')
            try:
                fp = fits.open(fname)
                fmt = 'fits'
            except IOError:
                pass

    # If we're here and fp is None, all attempt to open the file have
    # failed, so throw an IOError
    if fp is None:
        raise IOError(errno.ENOENT,
                      "unable to open file with base name "+fname)

    # Return the handle to the new file
    return fp, fname
