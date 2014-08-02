"""
Function to read a SLUG2 integrated_spec file.
"""

import numpy as np
from collections import namedtuple
import struct
from slug_open import slug_open

def read_integrated_spec(model_name, output_dir=None, asciionly=False,
                         binonly=False, verbose=False):
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
    verbose : bool
       If True, verbose output is printed as code runs

    Returns
    -------
    A namedtuple containing the following fields:
    wavelength : array
       wavelength, in Angstrom
    time : array
       times at which spectra are output, in yr
    L_lambda : array, shape (len(lambda), len(times))
       specific luminosity at each wavelength and each time, in erg/s/A
    """
    
    # Open file
    fp = slug_open(model_name+"_integrated_spec", output_dir=output_dir,
                   asciionly=asciionly, binonly=binonly)

    # Print status
    if verbose:
        print("Reading integrated spectra for model "+model_name)

    # Read ASCII or binary
    if fp.mode == 'r':

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

            # Split up the line
            data = entry.split()
            time.append(float(data[0]))
            wavelength.append(float(data[1]))
            L_lambda.append(float(data[2]))

        # Figure out how to decompose the data into arrays
        time = np.array(time)
        wavelength = np.array(wavelength)
        L_lambda = np.array(L_lambda)
        dummy, idx = np.unique(time, return_index=True)
        if len(idx) > 1:
            nl = idx[1]
            nt = len(time)/nl
            time = time[::nl]
            wavelength = wavelength[0:nl]
            L_lambda = np.reshape(L_lambda, (nt, nl))

    else:

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

        # Parse into arrays
        time = np.array(data_list[::nl+1])
        L_lambda = np.array(
            [data_list[(nl+1)*i+1:(nl+1)*(i+1)] 
             for i in range(nchunk)])

    # Close file
    fp.close()

    # Build the namedtuple to hold output
    out_type = namedtuple('integrated_spec',
                          ['wavelength', 'time', 'L_lambda'])
    out = out_type(wavelength, time, L_lambda)

    # Return
    return out
    
