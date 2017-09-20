"""
Function to read a SLUG2 integrated_cloudylines file.
"""

from collections import namedtuple
import numpy as np
import struct
from ..slug_open import slug_open

def read_integrated_cloudylines(model_name, output_dir=None, fmt=None,
                                verbose=False, read_info=None):
    """
    Function to read a SLUG2 integrated_cloudylines file.

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
       cloudy_linelabel : array, dtype='S4', shape (N_lines)
          labels for the lines, following cloudy's 4 character line label
          notation
       cloudy_linewl : array, shape (N_lines)
          rest wavelength for each line, in Angstrom
       cloudy_linelum : array, shape (N_lines, N_times, N_trials)
          luminosity of each line at each time for each trial, in erg/s
    """

    # Open file
    fp, fname = slug_open(model_name+"_integrated_cloudylines", 
                          output_dir=output_dir,
                          fmt=fmt)
    if read_info is not None:
        read_info['fname'] = fname

    # Print status
    if verbose:
        print("Reading integrated cloudy line luminosities for "
              "model "+model_name)
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
        label = []
        lum = []
        trial = []

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
            trial.append(trialptr)
            time.append(float(entry[0:11]))
            label.append(entry[21:25])
            wavelength.append(float(entry[28:39]))
            lum.append(float(entry[42:53]))

        # Convert to arrays
        trial = np.array(trial, dtype='uint')
        time = np.array(time)
        label = np.array(label)
        wavelength = np.array(wavelength)
        lum = np.array(lum)

        # Figure out the number of wavelengths by finding the first
        # time a wavelength repeats. Truncate the wavelength, label,
        # and time arrays appropriately.
        repeats = np.where(wavelength == wavelength[0])[0]
        if len(repeats > 1):
            nl = repeats[1]
            label = label[:nl]
            wavelength = wavelength[:nl]
            time = time[::nl]
        else:
            nl = len(wavelength)
            time = [time[0]]

        # Figure out how many trials there are and reshape the time
        # array appropriately
        ntrial = len(np.unique(trial))
        ntime = len(time)//ntrial
        if np.amin(time[:ntime] == time[ntime:2*ntime]):
            time = time[:ntime]

        # Reshape the line luminosity array
        lum = np.transpose(np.reshape(lum, (ntrial, ntime, nl)))

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Read number of lines
        nl = int(fp.readline())

        # Read line labels
        label = []
        for i in range(nl):
            label.append(fp.readline().split()[0])

        # Read the line wavelengths
        data = fp.read(struct.calcsize('d')*nl)
        wavelength = np.array(struct.unpack('d'*nl, data))

        # Now read the rest of the file and convert to doubles
        data = fp.read()
        nchunk = len(data) / \
                 (struct.calcsize('L')+(nl+1)*struct.calcsize('d'))
        data_list = struct.unpack(('L'+'d'*(nl+1))*nchunk, data)

        # Figure out how many times we have, and get unique times
        trial = np.array(data_list[::nl+2])
        time = np.array(data_list[1::nl+2])
        ntrial = len(np.unique(trial))
        ntime = len(time)//ntrial
        if np.amin(time[:ntime] == time[ntime:2*ntime]):
            time = time[:ntime]

        # Put line luminosity into array
        lum = np.zeros((nl, ntime, ntrial))
        ptr = 0
        for i in range(ntrial):
            for j in range(ntime):
                recptr = ptr*(nl+2)+2
                lum[:,j,i] = np.array(data_list[recptr:recptr+nl])
                ptr = ptr+1

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'

        # Get the line labels and wavelengths from the first HDU
        label = fp[1].data.field('Line_Label')
        wavelength = fp[1].data.field('Wavelength')

        # Get time, trial, line luminosities from the second HDU
        trial = fp[2].data.field('Trial')
        time = fp[2].data.field('Time')
        lum = fp[2].data.field('Line_Luminosity')

        # Re-arrange data into desired shape
        ntrial = len(np.unique(trial))
        ntime = len(time)//ntrial
        if ntime > 1:
            if np.amin(time[:ntime] == time[ntime:2*ntime]):
                time = time[:ntime]
        lum \
            = np.transpose(
                np.reshape(lum, (ntrial, ntime, len(wavelength))))

    # Close the file
    fp.close()

    # Build the namedtuple to hold output
    out_type = namedtuple('integrated_cloudyline',
                          ['time', 'cloudy_linelabel', 'cloudy_linewl', 
                           'cloudy_linelum'])
    out = out_type(time, label, wavelength, lum)

    # Return
    return out
