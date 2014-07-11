"""
Function to read a SLUG2 cluster_prop file.
"""

import numpy as np
from collections import namedtuple
import struct
import os
import os.path as osp

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

    # Did we get a specific directory in which to look? If not, try
    # current directory
    if output_dir is None:
        outdir = "."
    else:
        outdir = output_dir

    # See if we have a text file
    if not binonly:
        fname = osp.join(outdir, model_name+'_cluster_prop.txt')
        try:
            fp = open(fname, 'r')
        except IOError:
            fp = None
    else:
        fp = None

    # If that failed, look for a binary file
    if fp is None:
        if not asciionly:
            fname = osp.join(outdir, model_name+'_cluster_prop.bin')
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
                             model_name+'_cluster_prop.txt')
            try:
                fp = open(fname, 'r')
            except IOError:
                pass
        if not asciionly and fp is None:
            fname = osp.join(outdir, 
                             model_name+'_cluster_prop.bin')
            try:
                fp = open(fname, 'rb')
            except IOError:
                pass

    # If we're here and fp is None, all attempt to open the file have
    # failed, so throw an IOError
    if fp is None:
        raise IOError("unable to open cluster_prop file")

    # Print status
    if verbose:
        print("Reading file "+fname)

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

        # Suck file into memory
        data = fp.read()

        # Interpret data
        nentry = len(data)/struct.calcsize('LddddddQd')
        data_list = struct.unpack('LddddddQd'*nentry, data)

        # Stick data into correctly-named lists
        cluster_id = data_list[0::9]
        time = data_list[1::9]
        form_time = data_list[2::9]
        lifetime = data_list[3::9]
        target_mass = data_list[4::9]
        actual_mass = data_list[5::9]
        live_mass = data_list[6::9]
        num_star = data_list[7::9]
        max_star_mass = data_list[8::9]

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

