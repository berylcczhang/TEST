"""
Function to read a SLUG2 cluster_spec file.
"""

import numpy as np
from collections import namedtuple
import struct
from slug_open import slug_open

def read_cluster_spec(model_name, output_dir=None, asciionly=False,
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
    id : array, dtype uint
       unique ID of cluster
    time : array
       times at which cluster spectra are output, in yr
    L_lambda : array, shape (len(id), len(wavelength))
       specific luminosity of each cluster at each wavelength, in erg/s/A
    """

    # Open file
    fp = slug_open(model_name+"_cluster_spec", output_dir=output_dir,
                   asciionly=asciionly, binonly=binonly)

    # Print status
    if verbose:
        print("Reading "+model_name)

    # Prepare storage
    cluster_id = []
    time = []
    wavelength = []
    L_lambda = []

    # Read ASCII or binary
    if fp.mode == 'r':

        # ASCII mode

        # Burn the three header lines
        fp.readline()
        fp.readline()
        fp.readline()

        # Read first line and store cluster data
        entry = fp.readline()
        data = entry.split()
        cluster_id.append(long(data[0]))
        time.append(float(data[1]))
        wavelength.append(float(data[2]))
        L_lambda.append(float(data[3]))

        # Read the rest of the data for first cluster
        while True:
            entry = fp.readline()

            # Check for EOF
            if entry == '':
                break

            # Split up data
            data = entry.split()
            L_lambda.append(float(data[3]))
            id_tmp = long(data[0])

            # Stop when we find a different cluster
            if id_tmp != cluster_id[0]:
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
            data = entry.split()
            L_lambda.append(float(data[3]))
            ptr = ptr+1

            # When we get to the end of a chunk, push cluster ID and
            # time onto list, then reset pointer
            if ptr == nl:
                cluster_id.append(long(data[0]))
                time.append(float(data[1]))
                ptr = 0
                if verbose:
                    print("Read cluster {:d} at time {:e}".
                          format(cluster_id[-1], time[-1]))

    else:

        # Binary mode

        # First read number of wavelengths and wavelength table
        data = fp.read(struct.calcsize('L'))
        nl, = struct.unpack('L', data)
        data = fp.read(struct.calcsize('d')*nl)
        wavelength = np.array(struct.unpack('d'*nl, data))

        # Go through the rest of the file
        while True:

            # Read number of clusters and time in next block, checking
            # if we've hit eof
            data = fp.read(struct.calcsize('dL'))
            if len(data) < struct.calcsize('dL'):
                break
            t, ncluster = struct.unpack('dL', data)
            time.extend([t]*ncluster)

            # Read the next block of clusters
            data = fp.read(struct.calcsize('L')*ncluster + 
                           struct.calcsize('d')*ncluster*nl)
            data_list = struct.unpack(('L'+'d'*nl)*ncluster, data)

            # Pack clusters into data list
            cluster_id.extend(data_list[::nl+1])
            L_lambda.extend(
                [data_list[(nl+1)*i+1:(nl+1)*(i+1)] 
                 for i in range(ncluster)])

    # Close file
    fp.close()

    # Convert to arrays
    wavelength = np.array(wavelength)
    cluster_id = np.array(cluster_id, dtype='uint')
    time = np.array(time)
    L_lambda = np.array(L_lambda)
    L_lambda = np.reshape(L_lambda, (len(time), len(wavelength)))

    # Build namedtuple to hold output
    out_type = namedtuple('cluster_spec',
                          ['wavelength', 'id', 'time', 'L_lambda'])
    out = out_type(wavelength, cluster_id, time, L_lambda)

    # Return
    return out
