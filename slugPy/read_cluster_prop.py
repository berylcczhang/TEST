"""
Function to read a SLUG2 cluster_prop file.
"""

import numpy as np
from collections import namedtuple
import struct
from slug_open import slug_open

def read_cluster_prop(model_name, output_dir=None, fmt=None, 
                      verbose=False):
    """
    Function to read a SLUG2 integrated_prop file.

    Parameters
    ----------
    model_name : string
       The name of the model to be read
    output_dir : string
       The directory where the output is located; if set to None,
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
       number of living stars in cluster being treated stochastically
    max_star_mass : array
       mass of most massive living star in cluster
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_prop",
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading cluster properties for model "+model_name)

    # Prepare lists to hold data
    cluster_id = []
    trial = []
    time = []
    form_time = []
    lifetime = []
    target_mass = []
    actual_mass = []
    live_mass = []
    num_star = []
    max_star_mass = []

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode

        # Burn the three header lines
        fp.readline()
        fp.readline()
        fp.readline()

        # Read data
        trialptr = 0
        for entry in fp:
            if entry[:3] == '---':
                # Separator line
                trialptr = trialptr + 1
                continue
            data = entry.split()
            cluster_id.append(long(data[0]))
            trial.append(trialptr)
            time.append(float(data[1]))
            form_time.append(float(data[2]))
            lifetime.append(float(data[3]))
            target_mass.append(float(data[4]))
            actual_mass.append(float(data[5]))
            live_mass.append(float(data[6]))
            num_star.append(long(data[7]))
            max_star_mass.append(float(data[8]))

    elif fname.endswith('.bin'):

        # Binary mode

        # Set trial counter
        trialptr = 0

        # Go through file
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

            # Read the next block of clusters
            data = fp.read(struct.calcsize('LdddddQd')*ncluster)
            data_list = struct.unpack('LdddddQd'*ncluster, data)

            # If this time is not bigger than the last one was, this
            # is a new trial
            if len(time) > 0:
                if t <= time[-1]:
                    trialptr = trialptr + 1

            # Pack these clusters into the data list
            cluster_id.extend(data_list[0::8])
            time.extend([t]*ncluster)
            trial.extend([trialptr]*ncluster)
            form_time.extend(data_list[1::8])
            lifetime.extend(data_list[2::8])
            target_mass.extend(data_list[3::8])
            actual_mass.extend(data_list[4::8])
            live_mass.extend(data_list[5::8])
            num_star.extend(data_list[6::8])
            max_star_mass.extend(data_list[7::8])

    elif fname.endswith('.fits'):

        # FITS mode
        cluster_id = fp[1].data.field('UniqueID')
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')
        form_time = fp[1].data.field('FormTime')
        lifetime = fp[1].data.field('Lifetime')
        target_mass = fp[1].data.field('TargetMass')
        actual_mass = fp[1].data.field('BirthMass')
        live_mass = fp[1].data.field('LiveMass')
        num_star = fp[1].data.field('NumStar')
        max_star_mass = fp[1].data.field('MaxStarMass')

    # Close file
    fp.close()

    # Convert lists to arrays
    cluster_id = np.array(cluster_id, dtype='uint')
    trial = np.array(trial, dtype='uint')
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
                          ['id', 'trial', 'time', 'form_time', 
                           'lifetime', 'target_mass', 'actual_mass', 
                           'live_mass', 'num_star', 'max_star_mass'])
    out = out_type(cluster_id, trial, time, form_time, lifetime, 
                   target_mass, actual_mass, live_mass, num_star,
                   max_star_mass)

    # Return
    return out

