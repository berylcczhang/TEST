"""
Helper function to open slug output files.
"""

import os
import os.path as osp

def slug_open(filename, output_dir=None, asciionly=False, binonly=False):
    """
    Function to open a SLUG2 output file. By default this function
    searches both the current working directory and the output
    sub-directory of the SLUG_DIR directory, and it looks for both
    binary and text output.

    Parameters
    ----------
    filename : string
       Name of the file to open, without any extension
    output_dir : string
       The directory where the SLUG2 output is located; if set to None,
       the current directory is searched, followed by the SLUG_DIR
       directory if that environment variable is set
    asciionly : bool
       If True, only look for ASCII versions of outputs, ending in .txt
    binonly : bool
       If True, only look for binary versions of outputs, ending in .bin

    Returns
    -------
    fp : file
       A file object pointing the file that has been opened

    Raises
    ------
    IOError, if a file of the specified name cannot be found
    """

    # Did we get a specific directory in which to look? If not, try
    # current directory
    if output_dir is None:
        outdir = "."
    else:
        outdir = output_dir

    # See if we have a text file
    if not binonly:
        fname = osp.join(outdir, filename+'.txt')
        try:
            fp = open(fname, 'r')
        except IOError:
            fp = None
    else:
        fp = None

    # If that failed, look for a binary file
    if fp is None:
        if not asciionly:
            fname = osp.join(outdir, filename+'.bin')
            try:
                fp = open(fname, 'rb')
            except IOError:
                pass

    # If that failed, and we didn't get an explicit directory
    # specification, try looking in SLUG_DIR/output
    if (fp is None) and (output_dir is None) and \
       ('SLUG_DIR' in os.environ):
        outdir = osp.join(os.environ['SLUG_DIR'], 'output')
        if not binonly:
            fname = osp.join(outdir,
                             filename+'.txt')
            try:
                fp = open(fname, 'r')
            except IOError:
                pass
        if not asciionly and fp is None:
            fname = osp.join(outdir, 
                             filename+'.bin')
            try:
                fp = open(fname, 'rb')
            except IOError:
                pass

    # If we're here and fp is None, all attempt to open the file have
    # failed, so throw an IOError
    if fp is None:
        raise IOError("unable to open file with base name "+fname)

    # Return the handle to the new file
    return fp
