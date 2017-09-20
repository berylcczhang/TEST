"""
Function to read a SLUG2 integrated_cloudyphot file.
"""

import numpy as np
from collections import namedtuple
import struct
from ..slug_open import slug_open
from ..read_filter import read_filter

def read_integrated_cloudyphot(model_name, output_dir=None, fmt=None,
                               nofilterdata=False, photsystem=None,
                               verbose=False, read_info=None):
    """
    Function to read a SLUG2 integrated_cloudyphot file.

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

       time : array, shape (N_times) or shape (N_trials)
          Times at which data are output; shape is either N_times (if
          the run was done with fixed output times) or N_trials (if
          the run was done with random output times)
       cloudy_filter_names : list of string
          a list giving the name for each filter
       cloudy_filter_units : list of string
          a list giving the units for each filter
       cloudy_filter_wl_eff : list
          effective wavelength of each filter; this is set to None for the
          filters Lbol, QH0, QHe0, and QHe1; omitted if nofilterdata is
          True
       cloudy_filter_wl : list of arrays
          a list giving the wavelength table for each filter; this is
          None for the filters Lbol, QH0, QHe0, and QHe1; omitted if
          nofilterdata is True
       cloudy_filter_response : list of arrays
          a list giving the photon response function for each filter;
          this is None for the filters Lbol, QH0, QHe0, and QHe1; omitted
          if nofilterdata is True 
       cloudy_filter_beta : list
          powerlaw index beta for each filter; used to normalize the
          photometry
       cloudy_filter_wl_c : list
          pivot wavelength for each filter; used to normalize the photometry
       cloudy_phot_trans : array, shape (N_filter, N_times, N_trials)
          photometric value in each filter at each time in each trial for
          the transmitted light (i.e., the starlight remaining after it
          has passed through the HII region); units are as indicated in
          the units field
       cloudy_phot_emit : array, shape (N_filter, N_times, N_trials)
          photometric value in each filter at each time in each trial for
          the emitted light (i.e., the diffuse light emitted by the HII
          region); units are as indicated in the units field
       cloudy_phot_trans_emit : array, shape (N_filter, N_times, N_trials)
          photometric value in each filter at each time in each trial for
          the transmitted plus emitted light (i.e., the light coming
          directly from the stars after absorption by the HII region,
          plus the diffuse light emitted by the HII region); units are as
          indicated in the units field 

    Raises
       IOError, if no photometry file can be opened;
       ValueError, if photsystem is set to an unknown value
    """

    # Open file
    fp, fname = slug_open(model_name+"_integrated_cloudyphot", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading integrated cloudy photometry for model " +
              model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Read the list of filters
        line = fp.readline()
        filters = line.split()[1::3]
        nfilter = len(filters)

        # Burn a line
        line = fp.readline()

        # Read the list of units
        line = fp.readline()
        line = line.replace(')', '(').split('(') # split by ( and )
        units = []
        for l in line:
            if (not l.isspace()) and (len(l) > 0):
                units.append(l)
        units = units[1::3]   # Get rid of the units for time

        # Burn a line
        line = fp.readline()

        # Prepare holders for data
        trial = []
        time = []
        phot_trans = []
        phot_emit = []
        phot_trans_emit = []

        # Read through data
        trialptr = 0
        for line in fp:
            if line[:3] == '---':
                trialptr = trialptr + 1
                continue       # Skip separator lines
            linesplit = line.split()
            trial.append(trialptr)
            time.append(float(linesplit[0]))
            phot_trans.append(np.array(linesplit[1::3], dtype='float'))
            phot_emit.append(np.array(linesplit[2::3], dtype='float'))
            phot_trans_emit.append(np.array(linesplit[3::3], 
                                            dtype='float'))

        # Convert to arrays
        trial = np.array(trialptr)
        time = np.array(time)
        phot_trans = np.array(phot_trans)
        phot_emit = np.array(phot_emit)
        phot_trans_emit = np.array(phot_trans_emit)

        # Close
        fp.close()

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

        # Read rest of file, then close
        data = fp.read()
        fp.close()

        # Unpack the data
        chunkstr = 'L'+(nfilter+1)*'d'
        nchunk = len(data)/struct.calcsize(chunkstr)
        data_list = struct.unpack(nchunk*chunkstr, data)

        # Parse into arrays
        trial = np.array(data_list[::3*nfilter+2])
        time = np.array(data_list[1::3*nfilter+2])
        ntime = len(time)
        phot_trans = np.zeros((ntime, nfilter))
        phot_emit = np.zeros((ntime, nfilter))
        phot_trans_emit = np.zeros((ntime, nfilter))
        for i in range(ntime):
            ptr = (3*nfilter+2)*i
            phot_trans[i,:] = data_list[ptr+2:ptr+2+nfilter]
            phot_emit[i,:] = data_list[ptr+2+nfilter:ptr+2+2*nfilter]
            phot_trans_emit[i,:] = data_list[ptr+2+2*nfilter:ptr+2+3*nfilter]

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'

        # Get trial, time
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')

        # Get filter names and units
        filters = []
        units = []
        i = 3
        while 'TTYPE'+str(i) in fp[1].header.keys():
            filters.append(fp[1].header['TTYPE'+str(i)][:-12])
            units.append(fp[1].header['TUNIT'+str(i)])
            i = i+3

        # Get photometric data
        nfilter = len(filters)
        ntime = len(time)
        phot_trans = np.zeros((ntime, nfilter))
        phot_emit = np.zeros((ntime, nfilter))
        phot_trans_emit = np.zeros((ntime, nfilter))
        for i in range(len(filters)):
            phot_trans[:,i] = fp[1].data.field(filters[i]+'_Transmitted')
            phot_emit[:,i] = fp[1].data.field(filters[i]+'_Emitted')
            phot_trans_emit[:,i] \
                = fp[1].data.field(filters[i] + 
                                   '_Transmitted_plus_emitted')

    # Reshape time and photometry arrays
    ntrial = len(np.unique(trial))
    ntime = len(time)//ntrial
    if ntime > 1:
        if np.amin(time[:ntime] == time[ntime:2*ntime]):
            time = time[:ntime]
    phot_trans = np.transpose(np.reshape(phot_trans, 
                                         (ntrial, ntime, nfilter)))
    phot_emit = np.transpose(np.reshape(phot_emit, 
                                         (ntrial, ntime, nfilter)))
    phot_trans_emit = np.transpose(np.reshape(phot_trans_emit, 
                                         (ntrial, ntime, nfilter)))

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
            photometry_convert(photsystem, phot_trans, units)
            photometry_convert(photsystem, phot_emit, units)
            photometry_convert(photsystem, phot_trans_emit, units)
        else:
            photometry_convert(photsystem, phot_trans, units, wl_eff)
            photometry_convert(photsystem, phot_emit, units, wl_eff)
            photometry_convert(photsystem, phot_trans_emit, units, wl_eff)


    # Construct return object
    if nofilterdata:
        out_type = namedtuple('integrated_phot',
                              ['time', 'cloudy_filter_names', 
                               'cloudy_filter_units', 'cloudy_phot_trans',
                               'cloudy_phot_emit', 'cloudy_phot_trans_emit'])
        out = out_type(time, filters, units, phot_trans,
                       phot_emit, phot_trans_emit)
    else:
        out_type = namedtuple('integrated_phot',
                              ['time', 'cloudy_filter_names', 
                               'cloudy_filter_units',
                               'cloudy_filter_wl_eff', 
                               'cloudy_filter_wl', 
                               'cloudy_filter_response',
                               'cloudy_filter_beta',
                               'cloudy_filter_wl_c', 
                               'cloudy_phot_trans',
                               'cloudy_phot_emit', 
                               'cloudy_phot_trans_emit'])
        out = out_type(time, filters, units, wl_eff, wavelength, response,
                       beta, wl_c, phot_trans, phot_emit, phot_trans_emit)

    # Return
    return out
