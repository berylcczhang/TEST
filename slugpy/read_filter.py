"""
Function to read a filter response function for SLUG2
"""

import os
import os.path as osp
import numpy as np
from collections import namedtuple
import errno

def read_filter(filtername, filter_dir=None):
    """
    Function to read a filter or set of filters for SLUG2. By default
    this function searches the SLUG_DIR/lib/filter directory, followed
    by the current working directory. This can be overridden by the
    filter_dir keyword.

    Parameters
       filtername : string or iterable containing strings
          Name or names of filters to be read; for the special filters
          Lbol, QH0, QHe0, and QHe1, the return value will be None
       filter_dir : string
          Directory where the filter data files can be found

    Returns
       A namedtuple containing the following fields:

       wl_eff : float or array
          Central wavelength of the filter, defined by 
          wl_eff = exp(\int R ln lambda dln lambda / \int R dln lambda)
       wl : array or list of arrays
          Wavelength table for each filter, in Ang
       response : array or list of arrays
          Response function per photon for each filter
       beta : float or array
          Index beta for the filter
       wl_c : float or array
          Pivot wavelength for the filter; used when beta != 0 to
          normalize the photometry

    Raises
       IOError, if the filter data files cannot be opened, or if the
       requested filter cannot be found
    """

    # If filter list is not an iterable, make it an iterable of 1
    # element
    if type(filtername) is str:
        filter_names = [filtername]
    else:
        filter_names = filtername

    # If given an explicit directory, look in it
    if filter_dir is not None:

        fname1 = osp.join(filter_dir, 'FILTER_LIST')
        fname2 = osp.join(filter_dir, 'allfilters.dat')
        fp1 = open(fname1, 'r')
        fp2 = open(fname2, 'r')

    else:

        # No explicit directory given. If we have a SLUG_DIR
        # environment variable, look relative to it
        fp1 = None
        fp2 = None
        if 'SLUG_DIR' in os.environ:
            slugdir = os.environ['SLUG_DIR']
            fname1 = osp.join(slugdir, 'lib', 'filters',
                              'FILTER_LIST')
            fname2 = osp.join(slugdir, 'lib', 'filters',
                              'allfilters.dat')
            try:
                fp1 = open(fname1, 'r')
                fp2 = open(fname2, 'r')
            except IOError:
                # Close fp1 if we opened it and hit an error on fp2
                if fp1 is not None:
                    fp1.close()
                    fp1 = None

        # If we've failed to open, thus far, try the cwd
        if fp1 is None:
            fp1 = open('FILTER_LIST', 'r')
            fp2 = open('allfilters.dat', 'r')

    # If we're here, we've successfully opened both files. Now look
    # for the filter in the FILTER_LIST
    ctr = 0
    filteridx = [-2] * len(filter_names)
    filterbeta = [None] * len(filter_names)
    lambdac = [None] * len(filter_names)
    for line in fp1:

        # Get the filter name
        name = line.split()[1]

        # Compare name to list of filters we want; if there's a match,
        # record the index, beta, lambda_c
        if name in filter_names:
            filteridx[filter_names.index(name)] = ctr
            filterbeta[filter_names.index(name)] = float(line.split()[2])
            if line.split()[3] != '--':
                lambdac[filter_names.index(name)] = float(line.split()[3])

        # Increment counter
        ctr = ctr+1

    # Close the FILTER_LIST file
    fp1.close()

    # Set special filters to a flag value
    for i in range(len(filter_names)):
        if (filter_names[i] == 'QH0') or \
           (filter_names[i] == 'QHe0') or \
           (filter_names[i] == 'QHe1') or \
           (filter_names[i] == 'Lbol'):
            filteridx[i] = -1

    # Make sure we found a match for all filters; if not, throw an
    # error
    if -2 in filteridx:
        fp2.close()
        raise IOError(errno.EIO, "Unable to find filter " +
                      filter_names[filteridx.index(-2)] +
                      " in FILTER_LIST")

    # Burn the first line of allfilters.dat
    fp2.readline()

    # Create structure to hold results
    wavelength = [None] * len(filter_names)
    response = [None] * len(filter_names)

    # Now try to read the data from allfilters.dat
    ctr = 0
    wltmp = []
    rtmp = []
    for line in fp2:

        # Skip blank lines
        if len(line) == 0:
            continue

        # Is this the start of a new filter?
        if line.strip()[0] == '#':

            # Yes, this is a new filter. If we were recording at this
            # point, store what we've read in the result holder, then
            # reset the accumulators
            if ctr in filteridx:
                idx = filteridx.index(ctr)
                wavelength[idx] = np.array(wltmp)
                response[idx] = np.array(rtmp)
                wltmp = []
                rtmp = []

            # Increment the counter
            ctr = ctr+1

        else:

            # This is not the start of a new filter. If we're
            # recording this data, add it to the lists we're
            # building.
            if ctr in filteridx:
                wltmp.append(float(line.split()[0]))
                rtmp.append(float(line.split()[1]))

    # Close the file
    fp2.close()

    # If we were in the process of recording, store the last record
    if ctr in filteridx:
        idx = filteridx.index(ctr)
        wavelength[idx] = np.array(wltmp)
        response[idx] = np.array(rtmp)

    # Make sure we've successfully gotten data for all filters except
    # the special ones
    for i in range(len(wavelength)):
        if wavelength[i] is None and filteridx[i] != -1:
            raise IOError(errno.EIO, "Unable to find data for filter " +
                          filter_names[i] + " in allfilters.dat")

    # Compute the effective wavelength for each filter
    wl_eff = []
    for i in range(len(filter_names)):
        if filteridx[i] == -1:
            wl_eff.append(None)
        else:
            # Effective wavelength -- weighted by log lambda
            lnlambda = np.log(wavelength[i])
            dlnlambda = lnlambda[1:] - lnlambda[:-1]
            rlambda = 0.5*(lnlambda[1:]*response[i][1:] +
                           lnlambda[:-1]*response[i][:-1])
            rcen = 0.5 * (response[i][1:] + response[i][:-1])
            wl_eff.append(np.exp(np.sum(rlambda*dlnlambda) / 
                                 np.sum(rcen*dlnlambda)))
            # Commented out: weighted by lambda
            #dlambda = wavelength[i][1:] - wavelength[i][:-1]
            #rlambda = 0.5 * (wavelength[i][1:]*response[i][1:] + 
            #                 wavelength[i][:-1]*response[i][:-1])
            #rcen = 0.5 * (response[i][1:] + response[i][:-1])
            #wl_eff.append(np.sum(rlambda*dlambda) / np.sum(rcen*dlambda))

    # If we were given a scalar string as input, turn the output back
    # into scalars
    if type(filtername) is str:
        wavelength = wavelength[0]
        response = response[0]
        wl_eff = wl_eff[0]
        filterbeta = filterbeta[0]
        lambdac = lambdac[0]

    # Build the output object
    out_type = namedtuple('filter_data',
                          ['wl_eff', 'wl', 'response', 'beta', 'wl_c'])
    out = out_type(wl_eff, wavelength, response, filterbeta, lambdac)

    # Return
    return out
