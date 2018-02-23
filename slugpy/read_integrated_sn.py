"""
Function to read a SLUG2 integrated_sn file.
"""

import numpy as np
from collections import namedtuple
from collections import defaultdict
import struct
import re
from .slug_open import slug_open

def read_integrated_sn(model_name, output_dir=None, fmt=None, 
                       verbose=False, read_info=None):
    """
    Function to read a SLUG2 integrated_sn file.

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
    fp, fname = slug_open(model_name+"_integrated_sn", 
                          output_dir=output_dir,
                          fmt=fmt)

    # See if this file is a checkpoint file
    if len(re.findall('_chk\d\d\d\d', model_name)) != 0:
        checkpoint = True
    else:
        checkpoint = False
        
    # Print status and record
    if verbose:
        print("Reading integrated supernovae for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Prepare lists to hold data
    time = []
    tot_sn = []
    stoch_sn = []
    
    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # If this is a checkpoint file, skip the line stating how many
        # trials it contains; this line is not guaranteed to be
        # accurate, and is intended for the C++ code, not for us
        if checkpoint:
            fp.readline()

        # Read the first three header lines
        hdr = fp.readline()
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
            tot_sn.append(float(data[1]))
            stoch_sn.append(int(data[2]))

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # If this is a checkpoint, skip the bytes specifying how many
        # trials we have; this is inteded for the C++ code, not for
        # us, since we will determine that on our own
        if checkpoint:
            data = fp.read(struct.calcsize('i'))

        # Suck file into memory
        data = fp.read()

        # Interpret data
        datstr = 'LddL'
        nentry = len(data)//struct.calcsize(datstr)
        data_list = struct.unpack(datstr*nentry, data)

        # Stick data into correctly-named lists
        trial = data_list[0::4]
        time = data_list[1::4]
        tot_sn = data_list[2::4]
        stoch_sn = data_list[3::4]

    elif fname.endswith('fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')
        tot_sn = fp[1].data.field('TotSN')
        stoch_sn = fp[1].data.field('StochSN')

    # Close file
    fp.close()

    # Convert lists to arrays
    trial = np.array(trial)
    time = np.array(time)
    tot_sn = np.array(tot_sn)
    stoch_sn = np.array(stoch_sn, dtype='ulonglong')
    
    # Figure out if we have a number of trials with identical times,
    # indicating fixed output times, or if each trial has random times;
    # reshape time array appropriately
    ntrial = len(np.unique(trial))
    ntime = len(time)//ntrial
    if ntime != len(time):
        if np.amin(time[:ntime] == time[ntime:2*ntime]):
            time = time[:ntime]
            
    # Prune / reshape the output arrays
    tot_sn = np.transpose(tot_sn.reshape(ntrial, ntime))
    stoch_sn = np.transpose(stoch_sn.reshape(ntrial, ntime))
    
    # Build the namedtuple to hold output
    out_type = namedtuple('integrated_sn',
                          ['time', 'tot_sn', 'stoch_sn'])
    out = out_type(time, tot_sn, stoch_sn)

    # Return
    return out
