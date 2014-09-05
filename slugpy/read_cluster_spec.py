"""
Function to read a SLUG2 cluster_spec file.
"""

import numpy as np
from collections import namedtuple
import struct
from slug_open import slug_open

def read_cluster_spec(model_name, output_dir=None, fmt=None, 
                      verbose=False):
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

    Returns
    -------
    A namedtuple containing the following fields:
    id : array, dtype uint
       unique ID of cluster
    trial: array, dtype uint
       which trial was this cluster part of
    time : array
       times at which cluster spectra are output, in yr
    wl : array
       wavelength, in Angstrom
    spec : array, shape (N_cluster, N_wavelength)
       specific luminosity of each cluster at each wavelength, in erg/s/A

    Raises
    ------
    IOError, if no spectrum file can be opened
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_spec", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading cluster spectra for model "+model_name)

    # Prepare storage
    cluster_id = []
    time = []
    trial = []
    wavelength = []
    L_lambda = []

    # Read ASCII or binary
    if fname.endswith('.txt'):

        # ASCII mode

        # Burn the three header lines
        fp.readline()
        fp.readline()
        fp.readline()

        # Read first line and store cluster data
        trialptr = 0
        entry = fp.readline()
        data = entry.split()
        cluster_id.append(long(data[0]))
        time.append(float(data[1]))
        wavelength.append(float(data[2]))
        L_lambda.append(float(data[3]))
        trial.append(trialptr)

        # Read the rest of the data for first cluster
        while True:
            entry = fp.readline()

            # Check for EOF and separator lines
            if entry == '':
                break
            if entry[:3] == '---':
                trialptr = trialptr+1
                break

            # Split up data
            data = entry.split()
            L_lambda.append(float(data[3]))
            id_tmp = long(data[0])
            time_tmp = float(data[1])

            # Stop when we find a different cluster or a different time
            if id_tmp != cluster_id[0] or time_tmp != time[0]:
                break

            # Still the same cluster, so append to wavelength list
            wavelength.append(float(data[2]))

        # We have now read one full chunk, so we know how many
        # wavelength entries per cluster there are
        nl = len(wavelength)

        # Start of next chunk
        ptr = 1

        # Now read through rest of file
        while True:

            # Read a line
            entry = fp.readline()
            if entry == '':
                break
            if entry[:3] == '---':
                trialptr = trialptr+1
                continue
            data = entry.split()
            L_lambda.append(float(data[3]))
            ptr = ptr+1

            # When we get to the end of a chunk, push cluster ID,
            # time, trial number list, then reset pointer
            if ptr == nl:
                cluster_id.append(long(data[0]))
                time.append(float(data[1]))
                trial.append(trialptr)
                ptr = 0

    elif fname.endswith('.bin'):

        # Binary mode

        # First read number of wavelengths and wavelength table
        data = fp.read(struct.calcsize('L'))
        nl, = struct.unpack('L', data)
        data = fp.read(struct.calcsize('d')*nl)
        wavelength = np.array(struct.unpack('d'*nl, data))

        # Go through the rest of the file
        trialptr = 0
        while True:

            # Read number of clusters and time in next block, checking
            # if we've hit eof
            data = fp.read(struct.calcsize('dL'))
            if len(data) < struct.calcsize('dL'):
                break
            t, ncluster = struct.unpack('dL', data)

            # Skip if no clusters
            if ncluster == 0:
                continue

            # If this time is not bigger than the last one was, this
            # is a new trial
            if len(time) > 0:
                if t <= time[-1]:
                    trialptr = trialptr + 1

            # Add to time and trial arrays
            time.extend([t]*ncluster)
            trial.extend([trialptr]*ncluster)

            # Read the next block of clusters
            data = fp.read(struct.calcsize('L')*ncluster + 
                           struct.calcsize('d')*ncluster*nl)
            data_list = struct.unpack(('L'+'d'*nl)*ncluster, data)

            # Pack clusters into data list
            cluster_id.extend(data_list[::nl+1])
            L_lambda.extend(
                [data_list[(nl+1)*i+1:(nl+1)*(i+1)] 
                 for i in range(ncluster)])

    elif fname.endswith('.fits'):

        # FITS mode
        wavelength = fp[1].data.field('Wavelength')
        wavelength = wavelength.flatten()
        cluster_id = fp[2].data.field('UniqueID')
        trial = fp[2].data.field('Trial')
        time = fp[2].data.field('Time')
        L_lambda = fp[2].data.field('L_lambda')

    # Close file
    fp.close()

    # Convert to arrays
    wavelength = np.array(wavelength)
    cluster_id = np.array(cluster_id, dtype='uint')
    time = np.array(time)
    trial = np.array(trial, dtype='uint')
    L_lambda = np.array(L_lambda)
    L_lambda = np.reshape(L_lambda, (len(time), len(wavelength)))

    # Build namedtuple to hold output
    out_type = namedtuple('cluster_spec',
                          ['id', 'trial', 'time', 'wl', 'spec'])
    out = out_type(cluster_id, trial, time, wavelength, L_lambda)

    # Return
    return out
