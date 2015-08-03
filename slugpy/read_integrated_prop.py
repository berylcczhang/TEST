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
       fmt : 'txt' | 'ascii' | 'bin' | 'binary' | 'fits' | 'fits2'
          Format for the file to be read. If one of these is set, the
          function will only attempt to open ASCII-('txt' or 'ascii'), 
          binary ('bin' or 'binary'), or FITS ('fits' or 'fits2')
          formatted output, ending in .txt., .bin, or .fits,
          respectively. If set to None, the code will try to open
          ASCII files first, then if it fails try binary files, and if
          it fails again try FITS files.
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
       target_mass : array, shape (N_times, N_trials)
          Target stellar mass at each time
       actual_mass : array, shape (N_times, N_trials)
          Actual mass of stars created up to each time in each trial
       live_mass : array, shape (N_times, N_trials)
          Mass of currently-alive stars at each time in each trial
       stellar_mass : array
          mass of all stars, living and stellar remnants
       cluster_mass : array, shape (N_times, N_trials)
          Stellar mass in non-disrupted clusters at each time in each
          trial
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
    stellar_mass = []
    cluster_mass = []
    num_clusters = []
    num_dis_clusters = []
    num_fld_stars = []

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Read the first header line
        hdr = fp.readline()

        # See if we have the stellar mass field; this was added later,
        # so we check in order to maintain backwards compatibility
        hdrsplit = hdr.split()
        if 'StellarMass' in hdrsplit:
            has_st_mass = True
        else:
            has_st_mass = False

        # Burn two header lines
        fp.readline()
        fp.readline()

        # Read data
        trial = []
        trialptr = 0
        for entry in fp:
            if entry[:3] == '---':
                trialptr = trialptr+1
                continue       # Skip separator lines
            trial.append(trialptr)
            data = entry.split()
            time.append(float(data[0]))
            target_mass.append(float(data[1]))
            actual_mass.append(float(data[2]))
            live_mass.append(float(data[3]))
            if has_st_mass:
                stellar_mass.append(float(data[4]))
                cluster_mass.append(float(data[5]))
                num_clusters.append(int(data[6]))
                num_dis_clusters.append(int(data[7]))
                num_fld_stars.append(int(data[8]))
            else:
                stellar_mass.append(0.0)
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
        nentry = len(data)/struct.calcsize('LddddddQQQ')
        data_list = struct.unpack('LddddddQQQ'*nentry, data)

        # Stick data into correctly-named lists
        trial = data_list[0::10]
        time = data_list[1::10]
        target_mass = data_list[2::10]
        actual_mass = data_list[3::10]
        live_mass = data_list[4::10]
        stellar_mass = data_list[5::10]
        cluster_mass = data_list[6::10]
        num_clusters = data_list[7::10]
        num_dis_clusters = data_list[8::10]
        num_fld_stars = data_list[9::10]

    elif fname.endswith('fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')
        target_mass = fp[1].data.field('TargetMass')
        actual_mass = fp[1].data.field('ActualMass')
        live_mass = fp[1].data.field('LiveMass')
        try:
            stellar_mass = fp[1].data.field('StellarMass')
        except KeyError:
            stellar_mass = np.zeros(live_mass.shape)
        cluster_mass = fp[1].data.field('ClusterMass')
        num_clusters = fp[1].data.field('NumClusters')
        num_dis_clusters = fp[1].data.field('NumDisClust')
        num_fld_stars = fp[1].data.field('NumFldStar')

    # Close file
    fp.close()

    # Convert lists to arrays
    trial = np.array(trial)
    time = np.array(time)
    target_mass = np.array(target_mass)
    actual_mass = np.array(actual_mass)
    live_mass = np.array(live_mass)
    stellar_mass = np.array(stellar_mass)
    cluster_mass = np.array(cluster_mass)
    num_clusters = np.array(num_clusters, dtype='ulonglong')
    num_dis_clusters = np.array(num_dis_clusters, dtype='ulonglong')
    num_fld_stars = np.array(num_fld_stars, dtype='ulonglong')

    # Figure out if we have a number of trials with identical times,
    # indicating fixed output times, or if each trial has random times;
    # reshape time array appropriately
    ntrial = len(np.unique(trial))
    ntime = len(time)/ntrial
    if ntime != len(time):
        if np.amin(time[:ntime] == time[ntime:2*ntime]):
            time = time[:ntime]

    # Prune / reshape the output arrays
    target_mass = np.transpose(target_mass.reshape(ntrial, ntime))
    actual_mass = np.transpose(actual_mass.reshape(ntrial, ntime))
    live_mass = np.transpose(live_mass.reshape(ntrial, ntime))
    stellar_mass = np.transpose(stellar_mass.reshape(ntrial, ntime))
    cluster_mass = np.transpose(cluster_mass.reshape(ntrial, ntime))
    num_clusters = np.transpose(num_clusters.reshape(ntrial, ntime))
    num_dis_clusters = np.transpose(num_dis_clusters.reshape(ntrial, ntime))
    num_fld_stars = np.transpose(num_fld_stars.reshape(ntrial, ntime))

    # Build the namedtuple to hold output
    out_type = namedtuple('integrated_prop',
                          ['time', 'target_mass', 'actual_mass',
                           'live_mass', 'stellar_mass',
                           'cluster_mass', 'num_clusters',
                           'num_dis_clusters', 'num_fld_stars'])
    out = out_type(time, target_mass, actual_mass, live_mass,
                   stellar_mass, cluster_mass, num_clusters, 
                   num_dis_clusters, num_fld_stars)

    # Return
    return out

