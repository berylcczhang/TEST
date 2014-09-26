"""
Function to read a SLUG2 cluster_phot file.
"""

import numpy as np
from collections import namedtuple
from copy import deepcopy
import struct
from photometry_convert import photometry_convert
from read_filter import read_filter
from slug_open import slug_open

def read_cluster_phot(model_name, output_dir=None, fmt=None, 
                      nofilterdata=False, photsystem=None,
                      verbose=False, read_info=None):
    """
    Function to read a SLUG2 cluster_phot file.

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
       nofilterdata : bool
          If True, the routine does not attempt to read the filter
          response data from the standard location
       photsystem : None or string
          If photsystem is None, the data will be returned in the same
          photometric system in which they were read. Alternately, if it
          is a string, the data will be converted to the specified
          photometric system. Allowable values are 'L_nu', 'L_lambda',
          'AB', 'STMAG', and 'Vega', corresponding to the options defined
          in the SLUG code. If this is set and the conversion requested
          involves a conversion from a wavelength-based system to a
          frequency-based one, nofilterdata must be False so that the
          central wavelength of the photometric filters is available.
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
       filter_names : list of string
          a list giving the name for each filter
       filter_units : list of string
          a list giving the units for each filter
       filter_wl_eff : list
          effective wavelength of each filter; this is set to None for the
          filters Lbol, QH0, QHe0, and QHe1; omitted if nofilterdata is
          True
       filter_wl : list of arrays
          a list giving the wavelength table for each filter; this is
          None for the filters Lbol, QH0, QHe0, and QHe1; omitted if
          nofilterdata is True
       filter_response : list of arrays
          a list giving the photon response function for each filter;
          this is None for the filters Lbol, QH0, QHe0, and QHe1; omitted
          if nofilterdata is True 
       filter_beta : list
          powerlaw index beta for each filter; used to normalize the
          photometry
       filter_wl_c : list
          pivot wavelength for each filter; used to normalize the photometry
       phot : array, shape (N_cluster, N_filter)
          photometric value in each filter for each cluster; units are as
          indicated in the units field
       phot_ex : array, shape (N_filter, N_times, N_trials)
          same as phot, but after extinction has been applied (present
          only if SLUG was run with extinction enabled)
       
    Raises
       IOError, if no photometry file can be opened
       ValueError, if photsystem is set to an unknown values
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_phot", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading cluster photometry for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Prepare holders for data
    cluster_id = []
    time = []
    phot = []
    trial = []

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Read the list of filters
        line = fp.readline()
        filters = line.split()[2:]
        nfilter = len(filters)

        # Read the list of units
        line = fp.readline()
        line = line.replace(')', '(').split('(') # split by ( and )
        units = []
        for l in line:
            if (not l.isspace()) and (len(l) > 0):
                units.append(l)
        units = units[1:]    # Get rid of the units for time

        # See if we have extinction; this is indicated by there being
        # an even number of filters, and by the filters in the second
        # half of the list having the same names as those in the first
        # half, but with the extension "_ex"
        if nfilter % 2 == 0:
            extinct = True
            for i in range(nfilter/2):
                extinct = extinct and \
                          (filters[i]+'_ex' == filters[i+nfilter/2])
            phot_ex = []
        else:
            extinct = False

        # If we have extinction, reshape the filter and unit lists
        if extinct:
            nfilter = nfilter/2
            filters = filters[:nfilter]
            units = units[:nfilter]

        # Burn a line
        line = fp.readline()

        # Read through data
        trialptr = 0
        for line in fp:
            if line[:3] == '---':
                trialptr = trialptr+1
                continue
            linesplit = line.split()
            cluster_id.append(long(linesplit[0]))
            time.append(float(linesplit[1]))
            phot.append(linesplit[2:2+nfilter])
            if extinct:
                tmp_ex = []
                for i in range(nfilter):
                    substr = line[21*(2+nfilter+i):21*(2+nfilter+i+1)]
                    if substr.isspace():
                        tmp_ex.append(np.nan)
                    else:
                        tmp_ex.append(float(substr))
                phot_ex.append(np.array(tmp_ex))
            trial.append(trialptr)

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Read number of filters
        nfilter = int(fp.readline())

        # Read filter names and units
        filters = []
        units = []
        for i in range(nfilter):
            line = fp.readline()
            filters.append(line.split()[0])
            units.append(line.split()[1])

        # Read the bit that tells us if we're using extinction
        data = fp.read(struct.calcsize('b'))
        extinct = struct.unpack('b', data)[0] != 0
        if extinct:
            phot_ex = []

        # Go through the rest of the file
        trialptr = 0
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

            # If this time is not bigger than the last one was, this
            # is a new trial
            if len(time) > 0:
                if t <= time[-1]:
                    trialptr = trialptr + 1

            # Add to time and trial arrays
            time.extend([t]*ncluster)
            trial.extend([trialptr]*ncluster)

            # Read the next block of clusters
            if extinct:
                data = fp.read(struct.calcsize('L')*ncluster + 
                               struct.calcsize('d')*ncluster*nfilter*2)
                data_list = struct.unpack(('L'+'d'*nfilter*2)*ncluster, 
                                          data)
            else:
                data = fp.read(struct.calcsize('L')*ncluster + 
                               struct.calcsize('d')*ncluster*nfilter)
                data_list = struct.unpack(('L'+'d'*nfilter)*ncluster, data)

            # Pack clusters into data list
            if extinct:
                cluster_id.extend(data_list[::2*nfilter+1])
                phot.extend(
                    [data_list[(2*nfilter+1)*i+1:(2*nfilter+1)*i+1+nfilter] 
                     for i in range(ncluster)])
                phot_ex.extend(
                    [data_list[(2*nfilter+1)*i+1+nfilter:
                               (2*nfilter+1)*i+1+2*nfilter] 
                     for i in range(ncluster)])
            else:
                cluster_id.extend(data_list[::nfilter+1])
                phot.extend(
                    [data_list[(nfilter+1)*i+1:(nfilter+1)*(i+1)] 
                     for i in range(ncluster)])

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'

        # Get cluster ID, trial, time
        cluster_id = fp[1].data.field('UniqueID')
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')

        # Get filter names and units
        filters = []
        units = []
        i = 4
        while 'TTYPE'+str(i) in fp[1].header.keys():
            filters.append(fp[1].header['TTYPE'+str(i)])
            units.append(fp[1].header['TUNIT'+str(i)])
            i = i+1

        # See if we have extinction; this is indicated by there being
        # an even number of filters, and by the filters in the second
        # half of the list having the same names as those in the first
        # half, but with the extension "_ex"
        nfilter = len(filters)
        if nfilter % 2 == 0:
            extinct = True
            for i in range(nfilter/2):
                extinct = extinct and \
                          (filters[i]+'_ex' == filters[i+nfilter/2])
        else:
            extinct = False

        # If we have extinction, reshape the filter and unit lists
        if extinct:
            nfilter = nfilter/2
            filters = filters[:nfilter]
            units = units[:nfilter]

        # Get photometric data
        phot = np.zeros((len(time), len(filters)))
        if extinct:
            phot_ex = np.zeros((len(time), nfilter))
        for i in range(len(filters)):
            phot[:,i] = fp[1].data.field(filters[i])
            if extinct:
                phot_ex[:,i] = fp[1].data.field(filters[i]+"_ex")


    # Close file
    fp.close()

    # Convert to arrays
    cluster_id = np.array(cluster_id, dtype='uint')
    time = np.array(time, dtype='float')
    trial = np.array(trial, dtype='uint')
    phot = np.array(phot, dtype='float')
    phot = np.reshape(phot, (len(time), len(filters)))
    if extinct:
        phot_ex = np.array(phot_ex, dtype='float')
        phot_ex = np.reshape(phot_ex, (len(time), len(filters)))

    # Read filter data if requested
    if not nofilterdata:
        if verbose:
            print("Reading filter data")
        wl_eff, wavelength, response, beta, wl_c = read_filter(filters)

    # Do photometric system conversion if requested
    if photsystem is not None:
        if verbose:
            print("Converting photometric system")
        if nofilterdata:
            units_save = deepcopy(units)
            photometry_convert(photsystem, phot, units, filter_last=True)
            if extinct:
                photometry_convert(photsystem, phot_ex, units_save, 
                                   filter_names=filters)
        else:
            units_save = deepcopy(units)
            photometry_convert(photsystem, phot, units, wl_eff,
                               filter_last=True)
            if extinct:
                photometry_convert(photsystem, phot_ex, units_save, wl_eff, 
                                   filter_names=filters)

    # Construct return object
    if nofilterdata:
        if extinct:
            out_type = namedtuple('cluster_phot',
                                  ['id', 'trial', 'time', 'filter_names', 
                                   'filter_units', 'phot', 'phot_ex'])
            out = out_type(cluster_id, trial, time, filters, units, 
                           phot, phot_ex)
        else:
            out_type = namedtuple('cluster_phot',
                                  ['id', 'trial', 'time', 'filter_names', 
                                   'filter_units', 'phot'])
            out = out_type(cluster_id, trial, time, filters, units, phot)
    else:
        if extinct:
            out_type = namedtuple('cluster_phot',
                                  ['id', 'trial', 'time', 'filter_names', 
                                   'filter_units',
                                   'filter_wl_eff','filter_wl',
                                   'filter_response', 'filter_beta', 
                                   'filter_wl_c', 'phot', 'phot_ex'])
            out = out_type(cluster_id, trial, time, filters, units, wl_eff,
                           wavelength, response, beta, wl_c, phot, phot_ex)
        else:
            out_type = namedtuple('cluster_phot',
                                  ['id', 'trial', 'time', 'filter_names', 
                                   'filter_units',
                                   'filter_wl_eff','filter_wl',
                                   'filter_response', 'filter_beta', 
                                   'filter_wl_c', 'phot'])
            out = out_type(cluster_id, trial, time, filters, units, wl_eff,
                           wavelength, response, beta, wl_c, phot)

    # Return
    return out
