"""
Function to read a SLUG2 cluster_phot file.
"""

import numpy as np
from collections import namedtuple
import struct
from photometry_convert import photometry_convert
from read_filter import read_filter
from slug_open import slug_open

def read_cluster_phot(model_name, output_dir=None, asciionly=False,
                      binonly=False, nofilterdata=False, 
                      photsystem=None, verbose=False):
    """
    Function to read a SLUG2 integrated_phot file.

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

    Returns
    -------
    A namedtuple containing the following fields:
    filters : list of strings
       a list giving the name for each filter
    units : list of string
       a list giving the units for each filter
    central_wavelength : list
       central wavelength of each filter; this is set to None for the
       filters Lbol, QH0, QHe0, and QHe1; omitted if nofilterdata is
       True
    wavelength : list of arrays
       a list giving the wavelength table for each filter; this is
       None for the filters Lbol, QH0, QHe0, and QHe1; omitted if
       nofilterdata is True
    response : list of arrays
       a list giving the photon response function for each filter;
       this is None for the filters Lbol, QH0, QHe0, and QHe1; omitted
       if nofilterdata is True 
    id : array, dtype uint
       unique ID of cluster
    time : array
       times at which cluster photometry are output, in yr
    phot : array, shape (len(id), len(filters))
       photometric value in each filter for each cluster; units are as
       indicated in the units field
       
    Raises
    ------
    IOError, if no photometry file can be found
    ValueError, if photsystem is set to an unknown values
    """

    # Open file
    fp = slug_open(model_name+"_cluster_phot", output_dir=output_dir,
                   asciionly=asciionly, binonly=binonly)

    # Print status
    if verbose:
        print("Reading cluster photometry for model "+model_name)

    # Prepare holders for data
    cluster_id = []
    time = []
    phot = []

    # Read ASCII or binary
    if fp.mode == 'r':

        # ASCII mode

        # Read the list of filters
        line = fp.readline()
        filters = line.split()[2:]

        # Read the list of units
        line = fp.readline()
        line = line.replace(')', '(').split('(') # split by ( and )
        units = []
        for l in line:
            if (not l.isspace()) and (len(l) > 0):
                units.append(l)
        units = units[1:]    # Get rid of the units for time

        # Burn a line
        line = fp.readline()

        # Read through data
        for line in fp:
            linesplit = line.split()
            cluster_id.append(long(linesplit[0]))
            time.append(float(linesplit[1]))
            phot.append(linesplit[2:])

    else:

        # Binary mode

        # Read number of filters
        nfilter = int(fp.readline())

        # Read filter names and units
        filters = []
        units = []
        for i in range(nfilter):
            line = fp.readline()
            filters.append(line.split()[0])
            units.append(line.split()[1])

        # Go through the rest of the file
        while True:

            # Read number of clusters and time in next block, checking
            # if we've hit eof
            data = fp.read(struct.calcsize('dL'))
            if len(data) < struct.calcsize('dL'):
                break
            t, ncluster = struct.unpack('dL', data)
            time.extend([t]*ncluster)

            # Read the next block of clusters
            data = fp.read(struct.calcsize('L')*ncluster + 
                           struct.calcsize('d')*ncluster*nfilter)
            data_list = struct.unpack(('L'+'d'*nfilter)*ncluster, data)

            # Pack clusters into data list
            cluster_id.extend(data_list[::nfilter+1])
            phot.extend(
                [data_list[(nfilter+1)*i+1:(nfilter+1)*(i+1)] 
                 for i in range(ncluster)])

    # Close file
    fp.close()

    # Convert to arrays
    cluster_id = np.array(cluster_id, dtype='uint')
    time = np.array(time, dtype='float')
    phot = np.array(phot, dtype='float')
    phot = np.reshape(phot, (len(time), len(filters)))

    # Read filter data if requested
    if not nofilterdata:
        wl_cen, wavelength, response = read_filter(filters)

    # Do photometric system conversion if requested
    if photsystem is not None:
        if nofilterdata:
            photometry_convert(photsystem, phot, units)
        else:
            photometry_convert(photsystem, phot, units, wl_cen)

    # Construct return object
    if nofilterdata:
        out_type = namedtuple('cluster_phot',
                              ['filters', 'units', 'id', 'time', 'phot'])
        out = out_type(filters, units, clsuter_id, time, phot)
    else:
        out_type = namedtuple('integrated_phot',
                              ['filters', 'units',
                               'central_wavelength','wavelength',
                               'response', 'id', 'time', 'phot'])
        out = out_type(filters, units, wl_cen, wavelength, response,
                       cluster_id, time, phot)

    # Return
    return out
