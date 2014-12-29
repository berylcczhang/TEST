"""
Function to read a SLUG2 integrated_phot file.
"""

import numpy as np
from collections import namedtuple
from copy import deepcopy
import struct
from photometry_convert import photometry_convert
from read_filter import read_filter
from slug_open import slug_open

def read_integrated_phot(model_name, output_dir=None, fmt=None,
                         nofilterdata=False, photsystem=None,
                         verbose=False, read_info=None):
    """
    Function to read a SLUG2 integrated_phot file.

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
       phot : array, shape (N_filter, N_times, N_trials)
          photometric value in each filter at each time in each trial;
          units are as indicated in the units field
       phot_neb : array, shape (N_filter, N_times, N_trials)
          same as phot, but for the light after it has passed through
          the HII region (present only if SLUG was run with nebular
          emission enabled)
       phot_ex : array, shape (N_filter, N_times, N_trials)
          same as phot, but after extinction has been applied (present
          only if SLUG was run with extinction enabled)
       phot_neb_ex : array, shape (N_filter, N_times, N_trials)
          same as phot, but for the light after it has passed through
          the HII region and then had extinction applied (present only
          if SLUG was run with both nebular emission and extinction
          enabled)

    Raises
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
    if read_info is not None:
        read_info['fname'] = fname

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

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

        # Search for filters with names that end in _n, _ex, or _nex,
        # indicating that they include the effects of the nebula,
        # extinction, or both
        neb = []
        ex = []
        for i in range(nfilter):
            if len(filters[i]) > 2:
                if filters[i][-2:] == '_n':
                    neb.append(i)
            if len(filters[i]) > 3:
                if filters[i][-3:] == '_ex':
                    ext.append(i)
        if len(neb) > 0:
            nebular = True
        else:
            nebular = False
        if len(ex) > 0:
            extinct = True
        else:
            extinct = False

        # If we have nebular emission or extinction, reshape filter
        # and units lists
        nuniq = nfilter / ((1+nebular)*(1+extinct))
        nfilter = nuniq
        filters = filters[:nfilter]
        units = units[:nfilter]

        # Burn a line
        line = fp.readline()

        # Prepare holders for data
        trial = []
        time = []
        phot = []
        if nebular:
            phot_neb = []
        if extinct:
            phot_ex = []
            if nebular:
                phot_neb_ex = []

        # Read through data
        trialptr = 0
        for line in fp:
            if line[:3] == '---':
                trialptr = trialptr + 1
                continue       # Skip separator lines
            linesplit = line.split()
            trial.append(trialptr)
            time.append(float(linesplit[0]))
            phot.append(np.array(linesplit[1:nfilter+1],
                                 dtype='float'))
            if nebular:
                phot_neb.append(np.array(
                    linesplit[nfilter+1:2*nfilter+1],
                    dtype='float'))
            if extinct:
                tmp_ex = []
                for i in range(nfilter):
                    substr = line[21*(1+(1+nebular)*nfilter+i):
                                  21*(1+(1+nebular)*nfilter+i+1)]
                    if substr.isspace():
                        tmp_ex.append(np.nan)
                    else:
                        tmp_ex.append(float(substr))
                phot_ex.append(np.array(tmp_ex))
                if nebular:
                    tmp_neb_ex = []
                    for i in range(nfilter):
                        substr = line[21*(1+3*nfilter+i):
                                      21*(1+3*nfilter+i+1)]
                        if substr.isspace():
                            tmp_neb_ex.append(np.nan)
                        else:
                            tmp_neb_ex.append(float(substr))
                    phot_neb_ex.append(np.array(tmp_neb_ex))

        # Convert to arrays
        trial = np.array(trial)
        time = np.array(time)
        phot = np.array(phot)
        if nebular:
            phot_neb = np.array(phot_neb)
        if extinct:
            phot_ex = np.array(phot_ex)
            if nebular:
                phot_neb_ex = np.array(phot_neb_ex)

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

        # Read the bits that tells us if we're using nebular emission
        # and extinction
        data = fp.read(struct.calcsize('b'))
        nebular = struct.unpack('b', data)[0] != 0
        data = fp.read(struct.calcsize('b'))
        extinct = struct.unpack('b', data)[0] != 0

        # Read rest of file
        data = fp.read()

        # Unpack the data
        nftot = (1+nebular)*(1+extinct)*nfilter
        chunkstr = 'L'+(nftot+1)*'d'
        nchunk = len(data)/struct.calcsize(chunkstr)
        data_list = struct.unpack(nchunk*chunkstr, data)

        # Parse into arrays
        trial = np.array(data_list[::nftot+2], dtype=np.uint64)
        time = np.array(data_list[1::nftot+2])
        phot = np.zeros((nchunk, nfilter))
        if nebular:
            phot_neb = np.zeros((nchunk, nfilter))
        if extinct:
            phot_ex = np.zeros((nchunk, nfilter))
            if nebular:
                phot_neb_ex = np.zeros((nchunk, nfilter))
        for i in range(nchunk):
            phot[i,:] = data_list[(nftot+2)*i+2:
                                  (nftot+2)*i+2+nfilter]
            ptr = 1
            if nebular:
                phot_neb[i,:] = data_list[(nftot+2)*i+nfilter*ptr+2:
                                          (nftot+2)*i+nfilter*(ptr+1)+2]
                ptr = ptr+1
            if extinct:
                phot_ex[i,:] = data_list[(nftot+2)*i+nfilter*ptr+2:
                                         (nftot+2)*i+nfilter*(ptr+1)+2]
                ptr = ptr+1
                if nebular:
                    phot_neb_ex[i,:] \
                        = data_list[(nftot+2)*i+nfilter*ptr+2:
                                    (nftot+2)*i+nfilter*(ptr+1)+2]

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
            filters.append(fp[1].header['TTYPE'+str(i)])
            units.append(fp[1].header['TUNIT'+str(i)])
            i = i+1
        nfilter = len(filters)

        # Search for filters with names that end in _neb, _ex, or _neb_ex,
        # indicating that they include the effects of the nebula,
        # extinction, or both
        neb = []
        ex = []
        for i in range(nfilter):
            if len(filters[i]) > 4:
                if filters[i][-4:] == '_neb':
                    neb.append(i)
            if len(filters[i]) > 3:
                if filters[i][-3:] == '_ex':
                    ext.append(i)
        if len(neb) > 0:
            nebular = True
        else:
            nebular = False
        if len(ex) > 0:
            extinct = True
        else:
            extinct = False

        # If we have nebular emission or extinction, reshape filter
        # and units lists
        nuniq = nfilter / ((1+nebular)*(1+extinct))
        nfilter = nuniq
        filters = filters[:nfilter]
        units = units[:nfilter]

        # Get photometric data
        phot = np.zeros((len(time), nfilter))
        if nebular:
            phot_neb = np.zeros((len(time), nfilter))
        if extinct:
            phot_ex = np.zeros((len(time), nfilter))
            if nebular:
                phot_neb_ex = np.zeros((len(time), nfilter))
        for i in range(len(filters)):
            phot[:,i] = fp[1].data.field(filters[i])
            if nebular:
                phot_neb[:,i] = fp[1].data.field(filters[i]+"_neb")
            if extinct:
                phot_ex[:,i] = fp[1].data.field(filters[i]+"_ex")
                if nebular:
                    phot_neb_ex[:,i] = fp[1].data.field(filters[i]+"_neb_ex")

    # Close file
    fp.close()

    # Reshape time and photometry arrays
    ntrial = len(np.unique(trial))
    ntime = len(time)/ntrial
    if ntime != len(time):
        if np.amin(time[:ntime] == time[ntime:2*ntime]):
            time = time[:ntime]
    phot = np.transpose(np.reshape(phot, (ntrial, ntime, nfilter)))
    if nebular:
        phot_neb = np.transpose(np.reshape(phot_neb, (ntrial, ntime, nfilter)))
    if extinct:
        phot_ex = np.transpose(np.reshape(phot_ex, (ntrial, ntime, nfilter)))
        if nebular:
            phot_neb_ex = np.transpose(
                np.reshape(phot_neb_ex, (ntrial, ntime, nfilter)))

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
            photometry_convert(photsystem, phot, units, 
                               filter_names=filters)
            if nebular:
                photometry_convert(photsystem, phot_neb, units_save, 
                                   filter_names=filters)
            if extinct:
                photometry_convert(photsystem, phot_ex, units_save, 
                                   filter_names=filters)
                if nebular:
                    photometry_convert(photsystem, phot_neb_ex, units_save, 
                                       filter_names=filters)
        else:
            units_save = deepcopy(units)
            photometry_convert(photsystem, phot, units, wl_eff, 
                               filter_names=filters)
            units_out = deepcopy(units)
            units = deepcopy(units_save)
            if nebular:
                photometry_convert(photsystem, phot_neb, units, 
                                   wl_eff, filter_names=filters)
                units = deepcopy(units_save)
            if extinct:
                photometry_convert(photsystem, phot_ex, units, wl_eff, 
                                   filter_names=filters)
                units = deepcopy(units_save)
                if nebular:
                    photometry_convert(photsystem, phot_neb_ex, units, 
                                       wl_eff, filter_names=filters)
                    units = deepcopy(units_save)
            units = units_out

    # Construct return object
    fieldnames = ['time', 'filter_names', 'filter_units']
    fields = [time, filters, units]
    if not nofilterdata:
        fieldnames = fieldnames + ['filter_wl_eff', 'filter_wl',
                                   'filter_response', 'filter_beta',
                                   'filter_wl_c']
        fields = fields + [wl_eff, wavelength, response, beta, wl_c]
    fieldnames = fieldnames + ['phot']
    fields = fields + [phot]
    if nebular:
        fieldnames = fieldnames + ['phot_neb']
        fields = fields + [phot_neb]
    if extinct:
        fieldnames = fieldnames + ['phot_ex']
        fields = fields + [phot_ex]
        if nebular:
            fieldnames = fieldnames + ['phot_neb_ex']
            fields = fields + [phot_neb_ex]
    out_type = namedtuple('integrated_spec', fieldnames)
    out = out_type(*fields)

    # Return
    return out



        
