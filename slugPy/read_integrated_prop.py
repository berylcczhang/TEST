"""
Function to read a SLUG2 integrated_prop file.
"""

import numpy as np
from collections import namedtuple
import struct
import os
import os.path as osp

def read_integrated_prop(model_name, output_dir=None, asciionly=False,
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
    time : array, shape (N_trials, N_times)
       Output times, where N_times = number of output times, N_trials 
       = number of trials
    target_mass : array, shape (N_trials, N_times)
       Target stellar mass at each time
    actual_mass : array, shape (N_trials, N_times)
       Actual mass of stars created up to each time
    live_mass : array, shape (N_trials, N_times)
       Mass of currently-alive stars at each time
    cluster_mass : array, shape (N_trials, N_times)
       Mass of living stars in non-disrupted clusters at each time
    num_clusters : array, shape (N_trials, N_times)
       Number of non-disrupted clusters present at each time
    num_dis_clusters : array, shape (N_trials, N_times)
       Number of disrupted clusters present at each time
    num_fld_stars : array, shape (N_trials, N_times)
       Number of living field stars (excluding those in disrupted 
       clusters) present at each time
    """

    # Did we get a specific directory in which to look? If not, try
    # current directory
    if output_dir is None:
        outdir = "."
    else:
        outdir = output_dir

    # See if we have a text file
    if not binonly:
        fname = osp.join(outdir, model_name+'_integrated_prop.txt')
        try:
            fp = open(fname, 'r')
        except IOError:
            fp = None
    else:
        fp = None

    # If that failed, look for a binary file
    if fp is None:
        if not asciionly:
            fname = osp.join(outdir, model_name+'_integrated_prop.bin')
            try:
                fp = open(fname, 'rb')
            except IOError:
                pass

    # If that failed, and we didn't get an explicit directory
    # specification, try looking in SLUG_DIR/output
    if (fp is None) and (output_dir is None) and \
       ('SLUG_DIR' in os.environ):
        outdir = osp.join(os.environ['SLUG_DIR'], 'output')
        if not binonly:
            fname = osp.join(outdir,
                             model_name+'_integrated_prop.txt')
            try:
                fp = open(fname, 'r')
            except IOError:
                pass
        if not asciionly and fp is None:
            fname = osp.join(outdir, 
                             model_name+'_integrated_prop.bin')
            try:
                fp = open(fname, 'rb')
            except IOError:
                pass

    # If we're here and fp is None, all attempt to open the file have
    # failed, so throw an IOError
    if fp is None:
        raise IOError("unable to open integrated_prop file")

    # Print status
    if verbose:
        print("Reading file "+fname)

    # Prepare lists to hold data
    time = []
    target_mass = []
    actual_mass = []
    live_mass = []
    cluster_mass = []
    num_clusters = []
    num_dis_clusters = []
    num_fld_stars = []

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
            time.append(float(data[0]))
            target_mass.append(float(data[1]))
            actual_mass.append(float(data[2]))
            live_mass.append(float(data[3]))
            cluster_mass.append(float(data[4]))
            num_clusters.append(int(data[5]))
            num_dis_clusters.append(int(data[6]))
            num_fld_stars.append(int(data[7]))

    else:

        # Binary mode

        # Suck file into memory
        data = fp.read()

        # Interpret data
        nentry = len(data)/64
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

    # Close file
    fp.close()

    # Convert lists to arrays
    time = np.array(time)
    target_mass = np.array(target_mass)
    actual_mass = np.array(actual_mass)
    live_mass = np.array(live_mass)
    cluster_mass = np.array(cluster_mass)
    num_clusters = np.array(num_clusters, dtype='int')
    num_dis_clusters = np.array(num_dis_clusters, dtype='int')
    num_fld_stars = np.array(num_fld_stars, dtype='int')

    # Deduce number of times and trials by finding the first repeated
    # entry in the time data
    idx = np.where(time == time[0])
    ntrial = len(idx[0])
    ntime = len(time)/ntrial

    # Reshape the output arrays
    time = time.reshape((ntrial, ntime))
    target_mass = target_mass.reshape(ntrial, ntime)
    actual_mass = actual_mass.reshape(ntrial, ntime)
    live_mass = live_mass.reshape(ntrial, ntime)
    cluster_mass = cluster_mass.reshape(ntrial, ntime)
    num_clusters = num_clusters.reshape(ntrial, ntime)
    num_dis_clusters = num_dis_clusters.reshape(ntrial, ntime)
    num_fld_stars = num_fld_stars.reshape(ntrial, ntime)

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

