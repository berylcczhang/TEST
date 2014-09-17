"""
Function to read a SLUG2 integrated_prop file.
"""

import numpy as np
from collections import namedtuple
import struct
from slug_open import slug_open

def read_integrated_prop(model_name, output_dir=None, fmt=None, 
                         verbose=False, read_info=None):
    """
    Function to read a SLUG2 integrated_prop file.

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
          Times at which data are output
       target_mass : array, shape
          Target stellar mass at each time
       actual_mass : array, shape (N_times, N_trials)
          Actual mass of stars created up to each time in each trial
       live_mass : array, shape (N_times, N_trials)
          Mass of currently-alive stars at each time in each trial
       cluster_mass : array, shape (N_times, N_trials)
          Mass of living stars in non-disrupted clusters at each time in
          each trial
       num_clusters : array, shape (N_times, N_trials), dtype ulonglong
          Number of non-disrupted clusters present at each time in each
          trial
       num_dis_clusters : array, shape (N_times, N_trials), dtype ulonglong
          Number of disrupted clusters present at each time in each trial
       num_fld_stars : array, shape (N_times, N_trials), dtype ulonglong
          Number of living field stars (excluding those in disrupted 
          clusters and those being treated non-stochastically) present at
          each time in each trial
    """

    # Open file
    fp, fname = slug_open(model_name+"_integrated_prop", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status and record
    if verbose:
        print("Reading integrated properties for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Prepare lists to hold data
    time = []
    target_mass = []
    actual_mass = []
    live_mass = []
    cluster_mass = []
    num_clusters = []
    num_dis_clusters = []
    num_fld_stars = []

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Burn the three header lines
        fp.readline()
        fp.readline()
        fp.readline()

        # Read data
        for entry in fp:
            if entry[:3] == '---':
                continue       # Skip separator lines
            data = entry.split()
            time.append(float(data[0]))
            target_mass.append(float(data[1]))
            actual_mass.append(float(data[2]))
            live_mass.append(float(data[3]))
            cluster_mass.append(float(data[4]))
            num_clusters.append(int(data[5]))
            num_dis_clusters.append(int(data[6]))
            num_fld_stars.append(int(data[7]))

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Suck file into memory
        data = fp.read()

        # Interpret data
        nentry = len(data)/struct.calcsize('dddddQQQ')
        data_list = struct.unpack('dddddQQQ'*nentry, data)

        # Stick data into correctly-named lists
        time = data_list[0::8]
        target_mass = data_list[1::8]
        actual_mass = data_list[2::8]
        live_mass = data_list[3::8]
        cluster_mass = data_list[4::8]
        num_clusters = data_list[5::8]
        num_dis_clusters = data_list[6::8]
        num_fld_stars = data_list[7::8]

    elif fname.endswith('fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'
        time = fp[1].data.field('Time')
        target_mass = fp[1].data.field('TargetMass')
        actual_mass = fp[1].data.field('ActualMass')
        live_mass = fp[1].data.field('LiveMass')
        cluster_mass = fp[1].data.field('ClusterMass')
        num_clusters = fp[1].data.field('NumClusters')
        num_dis_clusters = fp[1].data.field('NumDisClust')
        num_fld_stars = fp[1].data.field('NumFldStar')

    # Close file
    fp.close()

    # Convert lists to arrays
    time = np.array(time)
    target_mass = np.array(target_mass)
    actual_mass = np.array(actual_mass)
    live_mass = np.array(live_mass)
    cluster_mass = np.array(cluster_mass)
    num_clusters = np.array(num_clusters, dtype='ulonglong')
    num_dis_clusters = np.array(num_dis_clusters, dtype='ulonglong')
    num_fld_stars = np.array(num_fld_stars, dtype='ulonglong')

    # Deduce number of times and trials by finding the first repeated
    # entry in the time data
    idx = np.where(time == time[0])
    ntrial = len(idx[0])
    ntime = len(time)/ntrial

    # Prune / reshape the output arrays
    time = time[0:ntime]
    target_mass = target_mass[0:ntime]
    actual_mass = np.transpose(actual_mass.reshape(ntrial, ntime))
    live_mass = np.transpose(live_mass.reshape(ntrial, ntime))
    cluster_mass = np.transpose(cluster_mass.reshape(ntrial, ntime))
    num_clusters = np.transpose(num_clusters.reshape(ntrial, ntime))
    num_dis_clusters = np.transpose(num_dis_clusters.reshape(ntrial, ntime))
    num_fld_stars = np.transpose(num_fld_stars.reshape(ntrial, ntime))

    # Build the namedtuple to hold output
    out_type = namedtuple('integrated_prop',
                          ['time', 'target_mass', 'actual_mass',
                           'live_mass', 'cluster_mass', 'num_clusters',
                           'num_dis_clusters', 'num_fld_stars'])
    out = out_type(time, target_mass, actual_mass, live_mass,
                   cluster_mass, num_clusters, num_dis_clusters,
                   num_fld_stars)

    # Return
    return out

