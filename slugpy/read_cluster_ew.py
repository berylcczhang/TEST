"""
Function to read a SLUG2 cluster_ew file.
"""

import numpy as np
from collections import namedtuple
from copy import deepcopy
import struct
import re
import errno
from .slug_open import slug_open

def read_cluster_ew(model_name, output_dir=None, fmt=None, 
                      verbose=False, read_info=None,
                      lines_only=False, read_lines=None,
                      ew_only=False):
    """
    Function to read a SLUG2 cluster_ew file.

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
          it fails again try FITS files. Only FITS files are supported
          for the Equivalent Width data file.
       verbose : bool
          If True, verbose output is printed as code runs
       read_info : dict
          On return, this dict will contain the keys 'fname' and
          'format', giving the name of the file read and the format it
          was in; 'format' will be one of 'ascii', 'binary', or 'fits'
       lines_only : bool
          If True, the code only reads the data on the lines, not
          any of the equivalent widths.
       read_lines : None | string | listlike containing strings
          If this is None, data on all lines is read. Otherwise only
          liness whose name(s) match the input line names are
          read.
       ew_only : bool
          If true, id, trial, time, and line information are not
          read, only the equivalent widths

    Returns
       A namedtuple, which can contain the following fields depending
       on the input options, and depending on which fields are present
       in the file being read:

       id : array, dtype uint
          unique ID of cluster
       trial: array, dtype uint
          which trial was this cluster part of
       time : array
          times at which cluster spectra are output, in yr
       line_names : list of string
          a list giving the name for each line
       line_units : list of string
          a list giving the units for the equivalent width of each line
       ew : array, shape (N_cluster, N_lines)
          equivalent width value of each line for each cluster; 
          units are as indicated in the units field
          
    Raises
       IOError, if no photometry file can be opened
       ValueError, if photsystem is set to an unknown value
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_ew", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading cluster equivalent widths for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # See if this file is a checkpoint file
    if len(re.findall('_chk\d\d\d\d', model_name)) != 0:
        checkpoint = True
    else:
        checkpoint = False

    # Prepare holders for data
    if not ew_only:
        cluster_id = []
        time = []
        trial = []
    ew = []

    # Read data
    if fname.endswith('.txt'):

        ########################################################
        # ASCII mode
        ########################################################

        raise NotImplementedError("EW only implemented for "
                                  "FITS files - ERROR - ABORTING")
        
    elif fname.endswith('.bin'):

        ########################################################
        # Binary mode
        ########################################################
        
        raise NotImplementedError("EW only implemented for "
                                  "FITS files - ERROR - ABORTING")
        

    elif fname.endswith('.fits'):

        ########################################################
        # FITS mode
        ########################################################

        if read_info is not None:
            read_info['format'] = 'fits'

        # Figure out which of the two fits formats we're using; if
        # there's > 2 HDUs, we're using fits2, otherwise we're using
        # fits1
        if len(fp) == 2:

            ########################################################
            # FITS1 format: all the data is stored in a single HDU
            ########################################################

            # Get line names and units
            lines = []
            units = []
            i = 4
            while 'TTYPE'+str(i) in fp[1].header.keys():
                lines.append(fp[1].header['TTYPE'+str(i)])
                units.append(fp[1].header['TUNIT'+str(i)])
                i = i+1
            nlines = len(lines)
            
            # If given a list of lines to read, make sure that we
            # haven't been given ones that are not available; if we have,
            # raise an error now
            if read_lines is not None:
                if hasattr(read_lines, '__iter__'):
                    for f in read_lines:
                        if not cl in lines:
                            raise IOError(errno.EIO,
                                          "requested line {:s} not available!".
                                          format(f))
                else:
                    if not read_lines in lines:
                        raise IOError(errno.EIO,
                                      "requested line {:s} not available!".
                                      format(read_lines))

 
            # If only reading lines, skip the rest
            if not lines_only:

                # If requested to read only certain lines, figure
                # out which ones are in our list to read
                if read_lines is not None:
                    lread = []
                    nl_final = 0
                    if hasattr(read_lines, '__iter__'):
                        for cl in lines:
                            if cl in read_lines:
                                lread.append(True)
                                nl_final = nl_final+1
                            else:
                                lread.append(False)
                    else:
                        for cl in lines:
                            if cl == read_lines:
                                lread.append(True)
                                nl_final = nl_final+1
                            else:
                                lread.append(False)
                else:
                    lread = [True] * len(lines)
                    nl_final = len(lines)

                # Get cluster ID, trial, time
                if not ew_only:
                    cluster_id = fp[1].data.field('UniqueID')
                    trial = fp[1].data.field('Trial')
                    time = fp[1].data.field('Time')

                # Get equivalent width data
                ew = np.zeros((fp[1].header['NAXIS2'], nl_final))

                ptr = 0
                for i in range(len(lines)):
                    if lread[i]:

                        ew[:,ptr] = fp[1].data.field(lines[i])
                        ptr = ptr+1

        else:

            ########################################################
            # FITS2 format: each line in its own HDU
            ########################################################

            # Get line names and units
            lines = [fp[i].header['TTYPE1'] for i in 
                       range(2,len(fp))]
            units = [fp[i].header['TUNIT1'] for i in 
                       range(2,len(fp))]                       
                       
            nlines = len(lines)

            # If given a list of lines to read, make sure that we
            # haven't been given ones that are not available; if we have,
            # raise an error now
            if read_lines is not None:
                if hasattr(read_lines, '__iter__'):
                    for cl in read_lines:
                        if not cl in lines:
                            raise IOError(errno.EIO,
                                          "requested line {:s} not available!".
                                          format(f))
                else:
                    if not read_lines in lines:
                        raise IOError(errno.EIO,
                                      "requested line {:s} not available!".
                                      format(read_lines))

          
            # If only reading lines, skip the rest
            if not lines_only:

               # If requested to read only certain lines, figure
                # out which ones are in our list to read
                if read_lines is not None:
                    lread = []
                    nl_final = 0
                    if hasattr(read_lines, '__iter__'):
                        for cl in lines:
                            if cl in read_lines:
                                lread.append(True)
                                nl_final = nl_final+1
                            else:
                                lread.append(False)
                    else:
                        for cl in lines:
                            if cl == read_lines:
                                lread.append(True)
                                nl_final = nl_final+1
                            else:
                                lread.append(False)
                else:
                    lread = [True] * len(lines)
                    nl_final = len(lines)

                # Get cluster ID, trial, time
                if not ew_only:
                    cluster_id = fp[1].data.field('UniqueID')
                    trial = fp[1].data.field('Trial')
                    time = fp[1].data.field('Time')

                # Get equivalent width

                ew = np.zeros((fp[1].header['NAXIS2'], nl_final))
                
                ptr1 = 0
                ptr2 = 2
                for i in range(len(lines)):
                    if lread[i]:
                        ew[:,ptr1] \
                            = fp[ptr2].data.field(lines[i])
                        ptr1 = ptr1 + 1
                    ptr2 = ptr2+1
                ptr1 = 0
                

    ########################################################
    # End of file reading block
    ########################################################

    # Close file
    fp.close()

    # If using only a subset of lines, truncate the line list now
    if read_lines is not None:
        lines_tmp = []
        units_tmp = []
        for i in range(len(lines)):
            if lread[i]:
                lines_tmp.append(lines[i])
                units_tmp.append(units[i])
        lines = lines_tmp
        units = units_tmp

    # Convert to arrays
    if not lines_only:
        if not ew_only:
            cluster_id = np.array(cluster_id, dtype='uint')
            time = np.array(time, dtype='float')
            trial = np.array(trial, dtype='uint')
        ew = np.array(ew, dtype='float')
        ew = np.reshape(ew, (ew.size/len(lines), 
                                 len(lines)))

    # Construct return object
    if lines_only or ew_only:
        fieldnames = ['line_names','line_units']
        fields = [lines,units]
    else:
        fieldnames = ['id', 'trial', 'time', 'line_names', 'line_units']
        fields = [cluster_id, trial, time, lines, units]


    fieldnames = fieldnames + ['ew']
    fields = fields + [ew]

    out_type = namedtuple('cluster_ew', fieldnames)
    out = out_type(*fields)

    # Return
    return out
