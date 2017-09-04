"""
Function to read a SLUG2 integrated_phot file.
"""

import numpy as np
from collections import namedtuple
from copy import deepcopy
import struct
import re
from .photometry_convert import photometry_convert
from .read_filter import read_filter
from .slug_open import slug_open

def read_integrated_phot(model_name, output_dir=None, fmt=None,
                         nofilterdata=False, photsystem=None,
                         verbose=False, read_info=None,
                         filters_only=False, read_filters=None,
                         read_nebular=None, read_extinct=None):
    """
    Function to read a SLUG2 integrated_phot file.

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
       filters_only : bool
          If True, the code only reads the data on the filters, not
          any of the actual photometry. If combined with nofilterdata,
          this can be used to return the list of available filters
          and nothing else.
       read_filters : None | string | listlike containing strings
          If this is None, data on all filters is read. Otherwise only
          filters whose name(s) match the input filter names ar
          read.
       read_nebular : None | bool
          If True, only data with the nebular contribution is read; if
          False, only data without it is read. Default behavior is to
          read all data.
       read_extinct : None | bool
          If True, only data with extinction applied is read; if
          False, only data without it is read. Default behavior is to
          read all data.

    Returns
       A namedtuple , which can contain the following fields depending
       on the input options, and depending on which fields are present
       in the file being read:

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
          None for the filters Lbol, QH0, QHe0, and QHe1
       filter_response : list of arrays
          a list giving the photon response function for each filter;
          this is None for the filters Lbol, QH0, QHe0, and QHe1
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
          the HII region
       phot_ex : array, shape (N_filter, N_times, N_trials)
          same as phot, but after extinction has been applied
       phot_neb_ex : array, shape (N_filter, N_times, N_trials)
          same as phot, but for the light after it has passed through
          the HII region and then had extinction applied

    Raises
       IOError, if no photometry file can be opened
       ValueError, if photsystem is set to an unknown value
    """

    # Open file
    fp, fname = slug_open(model_name+"_integrated_phot", 
                          output_dir=output_dir,
                          fmt=fmt)

    # See if this file is a checkpoint file
    if len(re.findall('_chk\d\d\d\d', model_name)) != 0:
        checkpoint = True
    else:
        checkpoint = False
        
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

        # If this is a checkpoint file, skip the line stating how many
        # trials it contains; this line is not guaranteed to be
        # accurate, and is intended for the C++ code, not for us
        if checkpoint:
            fp.readline()

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
                    ex.append(i)
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
        nuniq = nfilter // ((1+nebular)*(1+extinct))
        nfilter = nuniq
        filters = filters[:nfilter]
        units = units[:nfilter]

        # Skip the rest if reading filter data only
        if not filters_only:

            # If requested to read only certain filters, figure out
            # which ones are in our list to read
            if read_filters is not None:
                fread = []
                if hasattr(read_filters, '__iter__'):
                    for f in filters:
                        if f in read_filters:
                            fread.append(True)
                        else:
                            fread.append(False)
                else:
                    for f in filters:
                        if f == read_filters:
                            fread.append(True)
                        else:
                            fread.append(False)

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
            if nebular and extinct:
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
                if read_nebular is not True and read_extinct is not True:
                    if read_filters is None:
                        phot.append(np.array(linesplit[1:nfilter+1],
                                             dtype='float'))
                    else:
                        tmp = []
                        for i, j in enumerate(range(1,1+nfilter)):
                            if fread[i]:
                                tmp.append(linesplit[j])
                        phot.append(np.array(tmp, dtype='float'))
                if nebular and read_nebular is not False and \
                   read_extinct is not True:
                    if read_filters is None:
                        phot_neb.append(np.array(
                            linesplit[nfilter+1:2*nfilter+1],
                            dtype='float'))
                    else:
                        tmp = []
                        for i, j in enumerate(range(nfilter+1,1+2*nfilter)):
                            if fread[i]:
                                tmp.append(linesplit[j])
                        phot_neb.append(np.array(tmp, dtype='float'))
                if extinct and read_nebular is not True and \
                   read_extinct is not False:
                    tmp_ex = []
                    for i in range(nfilter):
                        substr = line[21*(1+(1+nebular)*nfilter+i):
                                      21*(1+(1+nebular)*nfilter+i+1)]
                        if read_filters is None:
                            if substr.isspace():
                                tmp_ex.append(np.nan)
                            else:
                                tmp_ex.append(float(substr))
                        else:
                            if fread[i]:
                                if substr.isspace():
                                    tmp_ex.append(np.nan)
                                else:
                                    tmp_ex.append(float(substr))
                    phot_ex.append(np.array(tmp_ex, dtype='float'))
                if nebular and extinct and read_nebular is not False \
                   and read_extinct is not False:
                    tmp_neb_ex = []
                    for i in range(nfilter):
                        substr = line[21*(1+3*nfilter+i):
                                      21*(1+3*nfilter+i+1)]
                        if read_filters is None:
                            if substr.isspace():
                                tmp_neb_ex.append(np.nan)
                            else:
                                tmp_neb_ex.append(float(substr))
                        else:
                            if fread[i]:
                                if substr.isspace():
                                    tmp_neb_ex.append(np.nan)
                                else:
                                    tmp_neb_ex.append(float(substr))
                    phot_neb_ex.append(np.array(tmp_neb_ex, dtype='float'))

            # Convert to arrays
            trial = np.array(trial)
            time = np.array(time)
            phot = np.array(phot)
            if nebular:
                phot_neb = np.array(phot_neb)
            if extinct:
                phot_ex = np.array(phot_ex)
            if nebular and extinct:
                phot_neb_ex = np.array(phot_neb_ex)

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # If this is a checkpoint, skip the bytes specifying how many
        # trials we have; this is inteded for the C++ code, not for
        # us, since we will determine that on our own
        if checkpoint:
            data = fp.read(struct.calcsize('i'))
        
        # Read number of filters
        nfilter = int(fp.readline())

        # Read filter names and units
        filters = []
        units = []
        for i in range(nfilter):
            line = fp.readline()
            filters.append(line.split()[0].decode("utf-8"))
            units.append(line.split()[1].decode("utf-8"))

        # Read the bits that tells us if we're using nebular emission
        # and extinction
        data = fp.read(struct.calcsize('b'))
        nebular = struct.unpack('b', data)[0] != 0
        data = fp.read(struct.calcsize('b'))
        extinct = struct.unpack('b', data)[0] != 0
        nftot = (1+nebular)*(1+extinct)*nfilter

        # If only reading the filters, skip the rest of this
        if not filters_only:

            # Handle things differently if we're reading
            # everything versus if we're reading just some filters
            if read_filters is None:

                # Read rest of file
                data = fp.read()

                # Unpack the data
                chunkstr = 'L'+(nftot+1)*'d'
                nchunk = len(data)//struct.calcsize(chunkstr)
                data_list = struct.unpack(nchunk*chunkstr, data)

                # Parse into arrays
                trial = np.array(data_list[::nftot+2], dtype=np.uint64)
                time = np.array(data_list[1::nftot+2])
                if read_nebular is not True and read_extinct is not True:
                    phot = np.zeros((nchunk, nfilter))
                if nebular and \
                   read_nebular is not False and \
                   read_extinct is not True:
                    phot_neb = np.zeros((nchunk, nfilter))
                if extinct and \
                   read_nebular is not True and \
                   read_extinct is not False:
                    phot_ex = np.zeros((nchunk, nfilter))
                if nebular and extinct and \
                   read_nebular is not False and \
                   read_extinct is not False:
                    phot_neb_ex = np.zeros((nchunk, nfilter))
                for i in range(nchunk):
                    if read_nebular is not True and \
                       read_extinct is not True:
                        phot[i,:] = data_list[(nftot+2)*i+2:
                                              (nftot+2)*i+2+nfilter]
                    ptr = 1
                    if nebular:
                        if read_nebular is not False and \
                           read_extinct is not True:
                            phot_neb[i,:] \
                                = data_list[(nftot+2)*i+nfilter*ptr+2:
                                            (nftot+2)*i+nfilter*(ptr+1)+2]
                        ptr = ptr+1
                    if extinct:
                        if read_nebular is not True and \
                           read_extinct is not False:
                            phot_ex[i,:] \
                                = data_list[(nftot+2)*i+nfilter*ptr+2:
                                            (nftot+2)*i+nfilter*(ptr+1)+2]
                        ptr = ptr+1
                    if nebular and extinct:
                        if read_nebular is not False and \
                           read_extinct is not False:
                            phot_neb_ex[i,:] \
                                = data_list[(nftot+2)*i+nfilter*ptr+2:
                                            (nftot+2)*i+nfilter*(ptr+1)+2]
                        ptr = ptr+1

            else:

                # Figure oput which filters to read
                fread = []
                if hasattr(read_filters, '__iter__'):
                    for f in filters:
                        if f in read_filters:
                            fread.append(True)
                        else:
                            fread.append(False)
                else:
                    for f in filters:
                        if f == read_filters:
                            fread.append(True)
                        else:
                            fread.append(False)

                # Set up lists to hold data
                trial = []
                time = []
                phot = []
                phot_neb = []
                phot_ex = []
                phot_neb_ex = []

                # Loop over times
                while True:

                    # Read trial and time
                    data = fp.read(struct.calcsize('Ld'))
                    if len(data) < struct.calcsize('Ld'):
                        break
                    data_list = struct.unpack('Ld', data)
                    trial.append(data_list[0])
                    time.append(data_list[1])

                    # Loop over filters
                    for i in range(nfilter):
                        if fread[i] and \
                           read_nebular is not True and \
                           read_extinct is not True:
                            data = fp.read(struct.calcsize('d'))
                            data_list = struct.unpack('d', data)
                            phot.extend(data_list)
                        else:
                            fp.seek(struct.calcsize('d'), 1)

                    # Repeat for nebular
                    if nebular:
                        for i in range(nfilter):
                            if fread[i] and \
                               read_nebular is not False and \
                               read_extinct is not True:
                                data = fp.read(struct.calcsize('d'))
                                data_list = struct.unpack('d', data)
                                phot_neb.extend(data_list)
                            else:
                                fp.seek(struct.calcsize('d'), 1)

                    # Repeat for extincted emission
                    if extinct:
                        for i in range(nfilter):
                            if fread[i] and \
                               read_nebular is not True and \
                               read_extinct is not False:
                                data = fp.read(struct.calcsize('d'))
                                data_list = struct.unpack('d', data)
                                phot_ex.extend(data_list)
                            else:
                                fp.seek(struct.calcsize('d'), 1)

                    # Repeat for nebular plus extincted emission
                    if nebular and extinct:
                        for i in range(nfilter):
                            if fread[i] and \
                               read_nebular is not False and \
                               read_extinct is not False:
                                data = fp.read(struct.calcsize('d'))
                                data_list = struct.unpack('d', data)
                                phot_neb_ex.extend(data_list)
                            else:
                                fp.seek(struct.calcsize('d'), 1)

                # Convert to arrays
                trial = np.array(trial, dtype=np.uint64)
                time = np.array(time, dtype='float')
                if read_nebular is not True and \
                   read_extinct is not True:
                    phot = np.array(phot, dtype='float')
                if nebular and \
                   read_nebular is not False and \
                   read_extinct is not True:
                    phot_neb = np.array(phot_neb, dtype='float')
                if extinct and \
                   read_nebular is not True and \
                   read_extinct is not False:
                    phot_ex = np.array(phot_ex, dtype='float')
                if nebular and extinct and \
                   read_nebular is not False and \
                   read_extinct is not False:
                    phot_neb_ex = np.array(phot_neb_ex, dtype='float')

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'

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
                    ex.append(i)
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
        nuniq = nfilter // ((1+nebular)*(1+extinct))
        nfilter = nuniq
        filters = filters[:nfilter]
        units = units[:nfilter]

        # If only reading filters, skip the rest
        if not filters_only:

            # If requested to read only certain filters, figure out
            # which ones are in our list to read
            if read_filters is not None:
                fread = []
                nf_final = 0
                if hasattr(read_filters, '__iter__'):
                    for f in filters:
                        if f in read_filters:
                            fread.append(True)
                            nf_final = nf_final+1
                        else:
                            fread.append(False)
                else:
                    for f in filters:
                        if f == read_filters:
                            fread.append(True)
                            nf_final = nf_final+1
                        else:
                            fread.append(False)
            else:
                fread = [True] * len(filters)
                nf_final = len(filters)

            # Get trial, time
            trial = fp[1].data.field('Trial')
            time = fp[1].data.field('Time')

            # Get photometric data
            if read_nebular is not True and \
               read_extinct is not True:
                phot = np.zeros((len(time), nf_final))
            if nebular and \
               read_nebular is not False and \
               read_extinct is not True:
                phot_neb = np.zeros((len(time), nf_final))
            if extinct and \
               read_nebular is not True and \
               read_extinct is not False:
                phot_ex = np.zeros((len(time), nf_final))
            if nebular and extinct and \
               read_nebular is not False and \
               read_extinct is not False:
                phot_neb_ex = np.zeros((len(time), nf_final))
            ptr = 0
            for i in range(len(filters)):
                if fread[i]:
                    if read_nebular is not True and \
                       read_extinct is not True:
                        phot[:,ptr] = fp[1].data.field(filters[i])
                    if nebular and \
                       read_nebular is not False and \
                       read_extinct is not True:
                        phot_neb[:,ptr] = fp[1].data.field(filters[i]+"_neb")
                    if extinct and \
                       read_nebular is not True and \
                       read_extinct is not False:
                        phot_ex[:,ptr] = fp[1].data.field(filters[i]+"_ex")
                    if nebular and extinct and \
                       read_nebular is not False and \
                       read_extinct is not False:
                        phot_neb_ex[:,ptr] \
                            = fp[1].data.field(filters[i]+"_neb_ex")
                    ptr = ptr+1

    # Close file
    fp.close()

    # If using only a subset of filters, truncate the filter list now
    if read_filters is not None:
        filters_tmp = []
        units_tmp = []
        for i in range(len(filters)):
            if fread[i]:
                filters_tmp.append(filters[i])
                units_tmp.append(units[i])
        filters = filters_tmp
        units = units_tmp
        nfilter = len(filters)

    # Reshape time and photometry arrays
    if not filters_only:
        ntrial = len(np.unique(trial))
        ntime = len(time)//ntrial
        if ntime != len(time):
            if np.amin(time[:ntime] == time[ntime:2*ntime]):
                time = time[:ntime]
        if read_nebular is not True and \
           read_extinct is not True:
            phot = np.transpose(np.reshape(phot, (ntrial, ntime, nfilter)))
        if nebular and \
           read_nebular is not False and \
           read_extinct is not True:
            phot_neb = np.transpose(np.reshape(phot_neb, 
                                               (ntrial, ntime, nfilter)))
        if extinct and \
           read_nebular is not True and \
           read_extinct is not False:
            phot_ex = np.transpose(np.reshape(phot_ex, 
                                              (ntrial, ntime, nfilter)))
        if nebular and extinct and \
           read_nebular is not False and \
           read_extinct is not False:
            phot_neb_ex = np.transpose(
                np.reshape(phot_neb_ex, (ntrial, ntime, nfilter)))

    # Read filter data if requested
    if not nofilterdata:
        if verbose:
            print("Reading filter data")
        wl_eff, wavelength, response, beta, wl_c = read_filter(filters)

    # Do photometric system conversion if requested
    if photsystem is not None and not filters_only:
        if verbose:
            print("Converting photometric system")
        if nofilterdata:
            units_save = deepcopy(units)
            if read_nebular is not True and \
               read_extinct is not True:
                photometry_convert(photsystem, phot, units, 
                                   filter_names=filters)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if nebular and \
               read_nebular is not False and \
               read_extinct is not True:
                photometry_convert(photsystem, phot_neb, units_save, 
                                   filter_names=filters)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if extinct and \
               read_nebular is not True and \
               read_extinct is not False:
                photometry_convert(photsystem, phot_ex, units_save, 
                                   filter_names=filters)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if nebular and extinct and \
               read_nebular is not False and \
               read_extinct is not False:
                photometry_convert(photsystem, phot_neb_ex, units_save, 
                                   filter_names=filters)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            units = units_out
        else:
            units_save = deepcopy(units)
            if read_nebular is not True and \
               read_extinct is not True:
                photometry_convert(photsystem, phot, units, wl_eff, 
                                   filter_names=filters)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if nebular and \
               read_nebular is not False and \
               read_extinct is not True:
                photometry_convert(photsystem, phot_neb, units, 
                                   wl_eff, filter_names=filters)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if extinct and \
               read_nebular is not True and \
               read_extinct is not False:
                photometry_convert(photsystem, phot_ex, units, wl_eff, 
                                   filter_names=filters)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if nebular and extinct and \
               read_nebular is not False and \
               read_extinct is not False:
                photometry_convert(photsystem, phot_neb_ex, units, 
                                   wl_eff, filter_names=filters)
                units = deepcopy(units_save)
                units_out = deepcopy(units)
            units = units_out

    # Construct return object
    if filters_only:
        fieldnames = ['filter_names', 'filter_units']
        fields = [filters, units]
    else:
        fieldnames = ['time', 'filter_names', 'filter_units']
        fields = [time, filters, units]
    if not nofilterdata:
        fieldnames = fieldnames + ['filter_wl_eff', 'filter_wl',
                                   'filter_response', 'filter_beta',
                                   'filter_wl_c']
        fields = fields + [wl_eff, wavelength, response, beta, wl_c]
    if not filters_only and read_nebular is not True and \
       read_extinct is not True:
        fieldnames = fieldnames + ['phot']
        fields = fields + [phot]
    if nebular and not filters_only and \
       read_nebular is not False and \
       read_extinct is not True:
        fieldnames = fieldnames + ['phot_neb']
        fields = fields + [phot_neb]
    if extinct and not filters_only and \
       read_nebular is not True and \
       read_extinct is not False:
        fieldnames = fieldnames + ['phot_ex']
        fields = fields + [phot_ex]
    if nebular and extinct and not filters_only and \
       read_nebular is not False and \
       read_extinct is not False:
        fieldnames = fieldnames + ['phot_neb_ex']
        fields = fields + [phot_neb_ex]
    out_type = namedtuple('integrated_phot', fieldnames)
    out = out_type(*fields)

    # Return
    return out



        
