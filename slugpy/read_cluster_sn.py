"""
Function to read a SLUG2 cluster_sn file.
"""

import numpy as np
from collections import namedtuple
import struct
import re
from .slug_open import slug_open

def read_cluster_sn(model_name, output_dir=None, fmt=None, 
                    verbose=False, read_info=None):
    """
    Function to read a SLUG2 cluster_sn file.

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

       id : array, dtype uint
          unique ID of cluster
       trial: array, dtype uint
          which trial was this cluster part of
       time : array
          times at which cluster spectra are output, in yr
       tot_sn : array, shape (N_cluster)
          total number of supernovae produced by each cluster up to
          the indicated time
       stoch_sn : array, shape (N_cluster)
          total number of supernovae produced by stars being treated
          stochastically in each cluster up to the indicated time

    Raises
       IOError, if no sn file can be opened
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_sn", 
                          output_dir=output_dir,
                          fmt=fmt)

    # See if this file is a checkpoint file
    if len(re.findall('_chk\d\d\d\d', model_name)) != 0:
        checkpoint = True
    else:
        checkpoint = False
        
    # Print status
    if verbose:
        print("Reading cluster supernovae for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Prepare storage
    cluster_id = []
    time = []
    trial = []
    tot_sn = []
    stoch_sn = []
    
    # Read ASCII or binary
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # If this is a checkpoint file, skip the line stating how many
        # trials it contains; this line is not guaranteed to be
        # accurate, and is intended for the C++ code, not for us
        if checkpoint:
            fp.readline()
        
        # Read the first header line
        hdr = fp.readline()

        # Burn the header lines
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
            cluster_id.append(int(data[0]))
            trial.append(trialptr)
            time.append(float(data[1]))
            tot_sn.append(float(data[2]))
            stoch_sn.append(int(data[3]))
            
    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # If this is a checkpoint, skip the bytes specifying how many
        # trials we have; this is inteded for the C++ code, not for
        # us, since we will determine that on our own
        if checkpoint:
            data = fp.read(struct.calcsize('i'))

        # Go through file
        while True:

            # Read number of clusters and time in next block, checking
            # if we've hit eof
            data = fp.read(struct.calcsize('LdL'))
            if len(data) < struct.calcsize('LdL'):
                break
            trialptr, t, ncluster = struct.unpack('LdL', data)

            # Skip if no clusters
            if ncluster == 0:
                continue

            # Read the next block of clusters
            datastr = 'LdL'
            data = fp.read(struct.calcsize(datastr)*ncluster)
            data_list = struct.unpack(datastr*ncluster, data)

            # Pack these clusters into the data list
            cluster_id.extend(data_list[0::3])
            time.extend([t]*ncluster)
            trial.extend([trialptr]*ncluster)
            tot_sn.extend(data_list[1::3])
            stoch_sn.extend(data_list[2::3])

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'
        cluster_id = fp[1].data.field('UniqueID')
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')
        tot_sn = fp[1].data.field('TotSN')
        stoch_sn = fp[1].data.field('StochSN')
        
    # Close file
    fp.close()

    # Convert lists to arrays
    cluster_id = np.array(cluster_id, dtype='uint')
    trial = np.array(trial, dtype='uint')
    time = np.array(time)
    tot_sn = np.array(tot_sn)
    stoch_sn = np.array(stoch_sn, dtype='uint')

    # Build the namedtuple to hold output
    out_list = ['id', 'trial', 'time', 'tot_sn', 
                'stoch_sn']
    out_dat = [cluster_id, trial, time, tot_sn, stoch_sn]
    out_type = namedtuple('cluster_int', out_list)
    out = out_type(*out_dat)
            
    # Return
    return out
