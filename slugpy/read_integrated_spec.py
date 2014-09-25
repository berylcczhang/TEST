"""
Function to read a SLUG2 integrated_spec file.
"""

import numpy as np
from collections import namedtuple
import struct
from slug_open import slug_open

def read_integrated_spec(model_name, output_dir=None, fmt=None, 
                         verbose=False, read_info=None):
    """
    Function to read a SLUG2 integrated_spec file.

    Parameters
       model_name : string
          The name of the model to be read
       output_dir : string
          The directory where the SLUG2 output is located; if set to None,
          the current directory is searched, followed by the SLUG_DIR
          directory if that environment variable is set
       fmt : string
          Format for the file to be read. Allowed values are 'ascii',
          'bin' or 'binary, and 'fits'. If one of these is set, the code
          will only attempt to open ASCII-, binary-, or FITS-formatted
          output, ending in .txt., .bin, or .fits, respectively. If set
          to None, the code will try to open ASCII files first, then if
          it fails try binary files, and if it fails again try FITS
          files.
       verbose : bool
          If True, verbose output is printed as code runs
       read_info : dict
          On return, this dict will contain the keys 'fname' and
          'format', giving the name of the file read and the format it
          was in; 'format' will be one of 'ascii', 'binary', or 'fits'

    Returns
       A namedtuple containing the following fields:

       time : array
          times at which spectra are output, in yr
       wl : array
          wavelength, in Angstrom
       spec : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity at each wavelength and each time for each
          trial, in erg/s/A
       wl_ex : array
          wavelength for the extincted spectrum, in Angstrom (present
          only if SLUG was run with extinction enabled)
       spec_ex : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity at each wavelength in wl_ex and each
          time for each trial, in erg/s/A (present only if SLUG was
          run with extinction enabled)
    """
    
    # Open file
    fp, fname = slug_open(model_name+"_integrated_spec", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading integrated spectra for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Prepare output holders
        wavelength = []
        time = []
        L_lambda = []

        # Read the first header line
        hdr = fp.readline()

        # See if we have extinction
        hdrsplit = hdr.split()
        if hdrsplit[-1] == 'L_lambda_ex':
            extinct = True
            wl_ex = []
            L_lambda_ex = []
        else:
            extinct = False

        # Burn the next two lines
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
            if len(data) > 3:
                wl_ex.append(float(data[1]))
                L_lambda_ex.append(float(data[3]))

        # Convert to arrays
        time = np.array(time)
        wavelength = np.array(wavelength)
        L_lambda = np.array(L_lambda)
        if extinct:
            wl_ex = np.array(wl_ex)
            L_lambda_ex = np.array(L_lambda_ex)

        # Figure out the number of wavelengths by finding the first
        # time a wavelength repeats. Truncate the wavelength and time
        # arrays appropriately.
        repeats = np.where(wavelength == wavelength[0])[0]
        if len(repeats) > 1:
            nl = repeats[1]
            wavelength = wavelength[:nl]
            time = time[::nl]
        else:
            nl = len(wavelength)
            time = [time[0]]
        if extinct:
            repeats = np.where(wl_ex == wl_ex[0])[0]
            if len(repeats) > 1:
                nl_ex = repeats[1]
                wl_ex = wl_ex[:nl_ex]

        # Figure out how many trials there are from how many times the
        # time array decreases instead of increasing. Truncate the
        # time array appropriately.
        ntrial = 1 + np.sum(time[1:] <= time[:-1])
        ntime = len(time)/ntrial
        time = time[:ntime]

        # Reshape the L_lambda array
        L_lambda = np.transpose(np.reshape(L_lambda, 
                                           (ntrial, ntime, nl)))
        if extinct:
            L_lambda_ex = np.transpose(
                np.reshape(L_lambda_ex, (ntrial, ntime, nl_ex)))

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Read a single character to see if extinction is included
        # in this file or not
        data = fp.read(struct.calcsize('b'))
        extinct = struct.unpack('b', data)[0] != 0

        # First read number of wavelengths and wavelength table
        data = fp.read(struct.calcsize('L'))
        nl, = struct.unpack('L', data)
        data = fp.read(struct.calcsize('d')*nl)
        wavelength = np.array(struct.unpack('d'*nl, data))
        if extinct:
            data = fp.read(struct.calcsize('L'))
            nl_ex, = struct.unpack('L', data)
            data = fp.read(struct.calcsize('d')*nl_ex)
            wl_ex = np.array(struct.unpack('d'*nl_ex, data))
        else:
            nl_ex = 0

        # Now read the rest of the file and convert to doubles
        data = fp.read()
        ndata = len(data)/struct.calcsize('d')
        nchunk = ndata/(nl+nl_ex+1)
        data_list = struct.unpack('d'*ndata, data)

        # Figure out how many times we have, and get unique times
        time = np.array(data_list[::nl+nl_ex+1])
        ntrial = 1 + np.sum(time[1:] <= time[:-1])
        ntime = len(time)/ntrial
        time = time[:ntime]

        # Put L_lambda into array
        L_lambda = np.zeros((nl, ntime, ntrial))
        if extinct:
            L_lambda_ex = np.zeros((nl_ex, ntime, ntrial))
        ptr = 0
        for i in range(ntrial):
            for j in range(ntime):
                ptr1 = ptr*(nl+nl_ex+1)+1
                L_lambda[:,j,i] \
                    = np.array(data_list[ptr1:ptr1+nl])
                if extinct:
                    L_lambda_ex[:,j,i] \
                        = np.array(data_list[ptr1+nl:ptr1+nl+nl_ex])
                ptr = ptr+1

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'

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

        # If we have extinction data, handle that too
        if 'Wavelength_ex' in fp[1].data.columns.names:
            extinct = True
            wl_ex = fp[1].data.field('Wavelength_ex')
            wl_ex = wl_ex.flatten()
            L_lambda_ex = fp[2].data.field('L_lambda_ex')
            L_lambda_ex \
                = np.transpose(
                    np.reshape(L_lambda_ex, (ntrial, ntime, len(wl_ex))))
        else:
            extinct = False

    # Close file
    fp.close()

    # Build the namedtuple to hold output
    if extinct:
        out_type = namedtuple('integrated_spec',
                              ['time', 'wl', 'spec', 'wl_ex', 'spec_ex'])
        out = out_type(time, wavelength, L_lambda, wl_ex, L_lambda_ex)
    else:
        out_type = namedtuple('integrated_spec',
                              ['time', 'wl', 'spec'])
        out = out_type(time, wavelength, L_lambda)

    # Return
    return out
    
