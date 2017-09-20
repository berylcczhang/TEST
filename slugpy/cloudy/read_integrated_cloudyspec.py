"""
Function to read a SLUG2 integrated_cloudyspec file.
"""

import numpy as np
from collections import namedtuple
import struct
from ..slug_open import slug_open

def read_integrated_cloudyspec(model_name, output_dir=None, fmt=None,
                               verbose=False, read_info=None):
    """
    Function to read a SLUG2 integrated_cloudyspec file.

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

       time : array, shape (N_times) or shape (N_trials)
          Times at which data are output; shape is either N_times (if
          the run was done with fixed output times) or N_trials (if
          the run was done with random output times)
       cloudy_wl : array
          wavelength, in Angstrom
       cloudy_inc : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity of the stellar radiation field at each
          wavelength and each time for each trial, in erg/s/A
       cloudy_trans : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity of the stellar radiation field after it has
          passed through the HII region, at each wavelength and each time
          for each trial, in erg/s/A
       cloudy_emit : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity of the radiation field emitted by the HII
          region, at each wavelength and each time for each trial, in
          erg/s/A
       cloudy_trans_emit : array, shape (N_wavelength, N_times, N_trials)
          the sum of emitted and transmitted; this is what would be seen
          by an observer looking at both the star cluster and its nebula
    """

    # Open file
    fp, fname = slug_open(model_name+"_integrated_cloudyspec", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading integrated cloudy spectra for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Prepare output holders
        trial = []
        wavelength = []
        time = []
        inc = []
        trans = []
        emit = []
        trans_emit = []

        # Burn the three header lines
        fp.readline()
        fp.readline()
        fp.readline()

        # Read data
        trialptr = 0
        for entry in fp:

            if entry[:3] == '---':
                trialptr = trialptr+1
                continue       # Skip separator lines

            # Split up the line
            data = entry.split()
            trial.append(trialptr)
            time.append(float(data[0]))
            wavelength.append(float(data[1]))
            inc.append(float(data[2]))
            trans.append(float(data[3]))
            emit.append(float(data[4]))
            trans_emit.append(float(data[5]))

        # Convert to arrays
        trial = np.array(trial)
        time = np.array(time)
        wavelength = np.array(wavelength)
        inc = np.array(inc)
        trans = np.array(trans)
        emit = np.array(emit)
        trans_emit = np.array(trans_emit)

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

        # Figure out how many trials there are and reshape the time
        # array appropriately
        ntrial = len(np.unique(trial))
        ntime = len(time)//ntrial
        if ntime > 1:
            if np.amin(time[:ntime] == time[ntime:2*ntime]):
                time = time[:ntime]

        # Reshape the spectral arrays
        inc = np.transpose(np.reshape(inc, (ntrial, ntime, nl)))
        trans = np.transpose(np.reshape(trans, (ntrial, ntime, nl)))
        emit = np.transpose(np.reshape(emit, (ntrial, ntime, nl)))
        trans_emit = np.transpose(np.reshape(trans_emit, (ntrial, ntime, nl)))

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # First read number of wavelengths and wavelength table
        data = fp.read(struct.calcsize('L'))
        nl, = struct.unpack('L', data)
        data = fp.read(struct.calcsize('d')*nl)
        wavelength = np.array(struct.unpack('d'*nl, data))

        # Now read the rest of the file and convert to correct type
        data = fp.read()
        nchunk = len(data) / \
                 (struct.calcsize('L')+(4*nl+1)*struct.calcsize('d'))
        data_list = struct.unpack(('L'+'d'*(4*nl+1))*nchunk, data)

        # Get time and trial arrays, and get number of times and trials
        trial = np.array(data_list[::4*nl+2], dtype='uint')
        time = np.array(data_list[1::4*nl+2])
        ntrial = len(np.unique(trial))
        ntime = len(time)//ntrial
        if ntime > 1:
            if np.amin(time[:ntime] == time[ntime:2*ntime]):
                time = time[:ntime]

        # Put spectra into arrays
        inc = np.zeros((nl, ntime, ntrial))
        trans = np.zeros((nl, ntime, ntrial))
        emit = np.zeros((nl, ntime, ntrial))
        trans_emit = np.zeros((nl, ntime, ntrial))
        ptr = 0
        for i in range(ntrial):
            for j in range(ntime):
                recptr = ptr*(4*nl+2)+2
                inc[:,j,i] = np.array(data_list[recptr:recptr+nl])
                trans[:,j,i] = np.array(data_list[recptr+nl:recptr+2*nl])
                emit[:,j,i] = np.array(data_list[recptr+2*nl:recptr+3*nl])
                trans_emit[:,j,i] = np.array(data_list[recptr+3*nl:recptr+4*nl])
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
        inc = fp[2].data.field('Incident_spectrum')
        trans = fp[2].data.field('Transmitted_spectrum')
        emit = fp[2].data.field('Emitted_spectrum')
        trans_emit = fp[2].data.field('Transmitted_plus_emitted_spectrum')

        # Re-arrange data into desired shape
        ntrial = len(np.unique(trial))
        ntime = len(time)//ntrial
        if ntime > 1:
            if np.amin(time[:ntime] == time[ntime:2*ntime]):
                time = time[:ntime]
        inc \
            = np.transpose(
                np.reshape(inc, (ntrial, ntime, len(wavelength))))
        trans \
            = np.transpose(
                np.reshape(trans, (ntrial, ntime, len(wavelength))))
        emit \
            = np.transpose(
                np.reshape(emit, (ntrial, ntime, len(wavelength))))
        trans_emit \
            = np.transpose(
                np.reshape(trans_emit, (ntrial, ntime, len(wavelength))))

    # Close the file
    fp.close()

    # Build the namedtuple to hold output
    out_type = namedtuple('integrated_cloudyspec',
                          ['time', 'cloudy_wl', 'cloudy_inc', 
                           'cloudy_trans', 'cloudy_emit', 
                           'cloudy_trans_emit'])
    out = out_type(time, wavelength, inc, trans, emit, trans_emit)

    # Return
    return out
    

