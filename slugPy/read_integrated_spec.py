"""
Function to read a SLUG2 integrated_spec file.
"""

import numpy as np
from collections import namedtuple
import struct
from slug_open import slug_open

def read_integrated_spec(model_name, output_dir=None, asciionly=False,
                         binonly=False, fitsonly=False, verbose=False):
    """
    Function to read a SLUG2 integrated_spec file.

    Parameters
    ----------
    model_name : string
       The name of the model to be read
    output_dir : string
       The directory where the SLUG2 output is located; if set to None,
       the current directory is searched, followed by the SLUG_DIR
       directory if that environment variable is set
    asciionly : bool
       If True, only look for ASCII versions of outputs, ending in .txt
    binonly : bool
       If True, only look for binary versions of outputs, ending in .bin
    fitsonly : bool
       If True, only look for FITS versions of outputs, ending in .fits
    verbose : bool
       If True, verbose output is printed as code runs

    Returns
    -------
    A namedtuple containing the following fields:
    time : array
       times at which spectra are output, in yr
    wl : array
       wavelength, in Angstrom
    spec : array, shape (N_wavelength, N_times, N_trials)
       specific luminosity at each wavelength and each time for each
       trial, in erg/s/A
    """
    
    # Open file
    fp, fname = slug_open(model_name+"_integrated_spec", 
                          output_dir=output_dir,
                          asciionly=asciionly, binonly=binonly,
                          fitsonly=fitsonly)

    # Print status
    if verbose:
        print("Reading integrated spectra for model "+model_name)

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode

        # Prepare output holders
        wavelength = []
        time = []
        L_lambda = []

        # Burn the three header lines
        fp.readline()
        fp.readline()
        fp.readline()

        # Read data
        for entry in fp:

            if entry[:3] == '---':
                continue       # Skip separator lines

            # Split up the line
            data = entry.split()
            time.append(float(data[0]))
            wavelength.append(float(data[1]))
            L_lambda.append(float(data[2]))

        # Convert to arrays
        time = np.array(time)
        wavelength = np.array(wavelength)
        L_lambda = np.array(L_lambda)

        # Figure out the number of wavelengths by finding the first
        # time a wavelength repeats. Truncate the wavelength and time
        # arrays appropriately.
        repeats = np.where(wavelength == wavelength[0])[0]
        if len(repeats > 1):
            nl = repeats[1]
            wavelength = wavelength[:nl]
            time = time[::nl]
        else:
            nl = len(wavelength)
            time = [time[0]]

        # Figure out how many trials there are from how many times the
        # time array decreases instead of increasing. Truncate the
        # time array appropriately.
        ntrial = 1 + np.sum(time[1:] <= time[:-1])
        ntime = len(time)/ntrial
        time = time[:ntime]

        # Reshape the L_lambda array
        L_lambda = np.transpose(np.reshape(L_lambda, (ntrial, ntime, nl)))

    elif fname.endswith('.bin'):

        # Binary mode

        # First read number of wavelengths and wavelength table
        data = fp.read(struct.calcsize('L'))
        nl, = struct.unpack('L', data)
        data = fp.read(struct.calcsize('d')*nl)
        wavelength = np.array(struct.unpack('d'*nl, data))

        # Now read the rest of the file and convert to doubles
        data = fp.read()
        ndata = len(data)/struct.calcsize('d')
        nchunk = ndata/(nl+1)
        data_list = struct.unpack('d'*ndata, data)

        # Figure out how many times we have, and get unique times
        time = np.array(data_list[::nl+1])
        ntrial = 1 + np.sum(time[1:] <= time[:-1])
        ntime = len(time)/ntrial
        time = time[:ntime]

        # Put L_lambda into array
        L_lambda = np.zeros((nl, ntime, ntrial))
        ptr = 0
        for i in range(ntrial):
            for j in range(ntime):
                L_lambda[:,j,i] \
                    = np.array(data_list[ptr*(nl+1)+1:(ptr+1)*(nl+1)])
                ptr = ptr+1

    elif fname.endswith('.fits'):

        # FITS mode

        # Read data
        wavelength = fp[1].data.field('Wavelength')
        wavelength = wavelength.flatten()
        trial = fp[2].data.field('Trial')
        time = fp[2].data.field('Time')
        L_lambda = fp[2].data.field('L_lambda')

        # Re-arrange data into desired shape
        ntrial = len(np.unique(trial))
        ntime = len(time)/ntrial
        time = time[:ntime]
        L_lambda \
            = np.transpose(
                np.reshape(L_lambda, (ntrial, ntime, len(wavelength))))

    # Close file
    fp.close()

    # Build the namedtuple to hold output
    out_type = namedtuple('integrated_spec',
                          ['time', 'wl', 'spec'])
    out = out_type(time, wavelength, L_lambda)

    # Return
    return out
    
