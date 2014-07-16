"""
Function to read a SLUG2 cluster_prop file.
"""

import numpy as np
from collections import namedtuple
import struct
from slug_open import slug_open

def read_cluster_prop(model_name, output_dir=None, asciionly=False,
                      binonly=False, verbose=False):
    """
    Function to read a SLUG2 integrated_prop file.

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
    id : array, dtype uint
       unique ID of cluster
    time : array
       time at which cluster's properties are being evaluated
    form_time : array
       time when cluster formed
    lifetime : array
       time at which cluster will disrupt
    target_mass : array
       target cluster mass
    actual_mass : array
       actual mass at formation
    live_mass : array
       mass of currently living stars
    num_star : array, dtype ulonglong
       number of living stars in cluster
    max_star_mass : array
       mass of most massive living star in cluster
    """

    # Open file
    fp = slug_open(model_name+"_cluster_prop", output_dir=output_dir,
                   asciionly=asciionly, binonly=binonly)

    # Print status
    if verbose:
        print("Reading "+model_name)

    # Prepare lists to hold data
    cluster_id = []
    time = []
    form_time = []
    lifetime = []
    target_mass = []
    actual_mass = []
    live_mass = []
    num_star = []
    max_star_mass = []

    # Read ASCII or binary
    if fp.mode == 'r':

        # ASCII mode

        # Burn the three header lines
        fp.readline()
        fp.readline()
        fp.readline()

        # Read data
        for entry in fp:
            data = entry.split()
            cluster_id.append(long(data[0]))
            time.append(float(data[1]))
            form_time.append(float(data[2]))
            lifetime.append(float(data[3]))
            target_mass.append(float(data[4]))
            actual_mass.append(float(data[5]))
            live_mass.append(float(data[6]))
            num_star.append(long(data[7]))
            max_star_mass.append(float(data[8]))

    else:

        # Binary mode

        # Go through file
        while True:

            # Read number of clusters and time in next block, checking
            # if we've hit eof
            data = fp.read(struct.calcsize('dL'))
            if len(data) < struct.calcsize('dL'):
                break
            t, ncluster = struct.unpack('dL', data)

            # Read the next block of clusters
            data = fp.read(struct.calcsize('LdddddQd')*ncluster)
            data_list = struct.unpack('LdddddQd'*ncluster, data)

            # Pack these clusters into the data list
            cluster_id.extend(data_list[0::8])
            time.extend([t]*ncluster)
            form_time.extend(data_list[1::8])
            lifetime.extend(data_list[2::8])
            target_mass.extend(data_list[3::8])
            actual_mass.extend(data_list[4::8])
            live_mass.extend(data_list[5::8])
            num_star.extend(data_list[6::8])
            max_star_mass.extend(data_list[7::8])

    # Close file
    fp.close()

    # Convert lists to arrays
    cluster_id = np.array(cluster_id, dtype='uint')
    time = np.array(time)
    form_time = np.array(form_time)
    lifetime = np.array(lifetime)
    target_mass = np.array(target_mass)
    actual_mass = np.array(actual_mass)
    live_mass = np.array(live_mass)
    num_star = np.array(num_star, dtype='ulonglong')
    max_star_mass = np.array(max_star_mass)

    # Build the namedtuple to hold output
    out_type = namedtuple('cluster_prop',
                          ['id', 'time', 'form_time', 'lifetime', 'target_mass',
                           'actual_mass', 'live_mass', 'num_star',
                           'max_star_mass'])
    out = out_type(cluster_id, time, form_time, lifetime, target_mass, 
                   actual_mass, live_mass, num_star, max_star_mass)

    # Return
    return out

