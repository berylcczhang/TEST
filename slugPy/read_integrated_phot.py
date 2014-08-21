"""
Function to read a SLUG2 integrated_phot file.
"""

import numpy as np
from collections import namedtuple
import struct
from photometry_convert import photometry_convert
from read_filter import read_filter
from slug_open import slug_open

def read_integrated_phot(model_name, output_dir=None, fmt=None,
                         nofilterdata=False, photsystem=None,
                         verbose=False):
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

    Returns
    -------
    A namedtuple containing the following fields:
    time : array
       times at which colors are output, in yr
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
    phot : array, shape (N_filter, N_times, N_trials)
       photometric value in each filter at each time in each trial;
       units are as indicated in the units field
       
    Raises
    ------
    IOError, if no photometry file can be opened
    ValueError, if photsystem is set to an unknown value
    """

    # Open file
    fp, fname = slug_open(model_name+"_integrated_phot", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading integrated photometry for model "+model_name)

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode

        # Read the list of filters
        line = fp.readline()
        filters = line.split()[1:]
        nfilter = len(filters)

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

        # Prepare holders for data
        time = []
        phot = []

        # Read through data
        for line in fp:
            if line[:3] == '---':
                continue       # Skip separator lines
            linesplit = line.split()
            time.append(float(linesplit[0]))
            phot.append(np.array(linesplit[1:], dtype='float'))

        # Convert to arrays
        time = np.array(time)
        phot = np.array(phot)

    elif fname.endswith('.bin'):

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

        # Read rest of file, then close
        data = fp.read()
        fp.close()

        # Unpack the data
        ndata = len(data)/struct.calcsize('d')
        ntime = ndata/(nfilter+1)
        data_list = struct.unpack('d'*ndata, data)

        # Parse into arrays
        time = np.array(data_list[::nfilter+1])
        phot = np.zeros((ntime, nfilter))
        for i in range(ntime):
            phot[i,:] = data_list[(nfilter+1)*i+1:(nfilter+1)*(i+1)]

    elif fname.endswith('.fits'):

        # FITS mode

        # Get trial, time
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')

        # Get filter names and units
        filters = []
        units = []
        i = 3
        while 'TTYPE'+str(i) in fp[1].header.keys():
            filters.append(fp[1].header['TTYPE'+str(i)])
            units.append(fp[1].header['TUNIT'+str(i)])
            i = i+1

        # Get photometric data
        nfilter = len(filters)
        phot = np.zeros((len(time), nfilter))
        for i in range(len(filters)):
            phot[:,i] = fp[1].data.field(filters[i])

    # Close file
    fp.close()

    # Reshape time and photometry arrays
    ntrial = 1 + np.sum(time[1:] <= time[:-1])
    ntime = len(time)/ntrial
    time = time[:ntime]
    phot = np.transpose(np.reshape(phot, (ntrial, ntime, nfilter)))

    # Read filter data if requested
    if not nofilterdata:
        if verbose:
            print("Reading filter data")
        wl_eff, wavelength, response = read_filter(filters)

    # Do photometric system conversion if requested
    if photsystem is not None:
        if verbose:
            print("Converting photometric system")
        if nofilterdata:
            photometry_convert(photsystem, phot, units)
        else:
            photometry_convert(photsystem, phot, units, wl_eff)

    # Construct return object
    if nofilterdata:
        out_type = namedtuple('integrated_phot',
                              ['time', 'filter_names', 
                               'filter_units', 'phot'])
        out = out_type(time, filters, units, phot)
    else:
        out_type = namedtuple('integrated_phot',
                              ['time', 'filter_names', 'filter_units',
                               'filter_wl_eff', 'filter_wl', 
                               'filter_response', 'phot'])
        out = out_type(time, filters, units, wl_eff, wavelength, response,
                       phot)

    # Return
    return out



        
