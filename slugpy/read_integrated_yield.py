"""
Function to read a SLUG2 integrated_yield file.
"""

import numpy as np
from collections import namedtuple
import struct
from slug_open import slug_open

def read_integrated_yield(model_name, output_dir=None, fmt=None,
                          read_info=None, verbose=False):

    """
    Function to read a SLUG2 integrated_yield file.

    Parameters
       model_name : string
          The name of the model to be read
       output_dir : string
          The directory where the SLUG2 output is located; if set to None,
          the current directory is searched, followed by the SLUG_DIR
          directory if that environment variable is set
       fmt : 'txt' | 'ascii' | 'bin' | 'binary' | 'fits' | 'fits2'
          Format for the file to be read. If one of these is set, the
          function will only attempt to open ASCII-('txt' or 'ascii'), 
          binary ('bin' or 'binary'), or FITS ('fits' or 'fits2')
          formatted output, ending in .txt., .bin, or .fits,
          respectively. If set to None, the code will try to open
          ASCII files first, then if it fails try binary files, and if
          it fails again try FITS files.
       verbose : bool
          If True, verbose output is printed as code runs
       read_info : dict
          On return, this dict will contain the keys 'fname' and
          'format', giving the name of the file read and the format it
          was in; 'format' will be one of 'ascii', 'binary', or 'fits'

    Returns
       A namedtuple containing the following fields:

       name : array of strings, shape (N_iso)
          Atomic symbols of the isotopes included in the yield table
       Z : array of int, shape (N_iso)
          Atomic numbers of the isotopes included in the yield table
       A : array of int, shape (N_iso)
          Atomic mass number of the isotopes included in the yield table
       time : array, shape (N_times) or shape (N_trials)
          Times at which data are output; shape is either N_times (if
          the run was done with fixed output times) or N_trials (if
          the run was done with random output times)
       yld : array, shape (N_times, N_iso) or (N_trials, N_iso)
          Yield of each isotope, defined as the instantaneous amount
          produced up to that time; for unstable isotopes, this
          includes the effects of decay since production
    """

    # Open file
    fp, fname = slug_open(model_name+"_integrated_yield", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status and record
    if verbose:
        print("Reading integrated yield for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Burn 3 header lines
        hdr = fp.readline()
        hdr = fp.readline()
        hdr = fp.readline()

        # Read data
        time = []
        isotope_name = []
        isotope_Z = []
        isotope_A = []
        yld = []
        yldtmp = []
        first_trial = True
        for entry in fp:
            # See if this is a new trial
            if entry[:3] == '---':
                yld.append(yldtmp)
                yldtmp = []
                first_trial = False
                continue       # Skip separator lines
            # Read data
            data = entry.split()
            time.append(float(data[0]))
            if first_trial:
                isotope_name.append(data[1].title())
                isotope_Z.append(float(data[2]))
                isotope_A.append(float(data[3]))
            yldtmp.append(float(data[4]))

        # Append last set of yields
        yld.append(yldtmp)

        # Truncate repeats in the isotope list that correspond to the
        # same set of isotopes at different times, and get unique times
        if len(isotope_name) > 1:
            isotopes = zip(isotope_name, isotope_Z, isotope_A)
            if isotopes[0] in isotopes[1:]:
                niso = isotopes[1:].index(isotopes[0])+1
                isotope_name = isotope_name[:niso]
                isotope_Z = isotope_Z[:niso]
                isotope_A = isotope_A[:niso]
                time = time[::niso]

        # Convert to arrays
        isotope_name = np.array(isotope_name)
        isotope_Z = np.array(isotope_Z, dtype=int)
        isotope_A = np.array(isotope_A, dtype=int)
        time = np.array(time)

        # If times are repeats, truncate
        if time.size > len(yld):
            ntime = time.size/len(yld)
            if np.amin(time[:ntime] == time[ntime:2*ntime]):
                time = time[:ntime]

        # Make yield array
        yld = np.array(yld)
        yld = yld.reshape((yld.shape[0], time.size,
                           isotope_name.size))

    # Close file
    fp.close()

    # Build output holder
    fieldnames = ['name', 'Z', 'A', 'time', 'yld']
    fields = [ isotope_name, isotope_Z, isotope_A, time, yld]
    out_type = namedtuple('integrated_yield', fieldnames)
    out = out_type(*fields)

    # Return
    return out


