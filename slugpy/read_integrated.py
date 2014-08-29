"""
Function to read all integrated data for a SLUG2 run.
"""

from collections import namedtuple
from read_integrated_prop import read_integrated_prop
from read_integrated_phot import read_integrated_phot
from read_integrated_spec import read_integrated_spec
from cloudy.read_integrated_cloudyspec import read_integrated_cloudyspec

def read_integrated(model_name, output_dir=None, fmt=None,
                    nofilterdata=False, photsystem=None, 
                    verbose=False, read_info=None):
    """
    Function to read all integrated data for a SLUG2 run.

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
    read_info : dict
       On return, this dict will contain the keys 'prop_name',
       'phot_name', 'spec_name', 'cloudyspec_name', and 'format',
       giving the names of the files read and the format they were in;
       'format' will be one of 'ascii', 'binary', or 'fits'

    Returns
    -------
    A namedtuple containing the following fields:

    (Always present)
    time: array
       Times at which data are output

    (Only present if an integrated_prop file is found)
    target_mass : array, shape (N_times)
       Target stellar mass at each time
    actual_mass : array, shape (N_times, N_trials)
       Actual mass of stars created up to each time in each trial
    live_mass : array, shape (N_times, N_trials)
       Mass of currently-alive stars at each time in each trial
    cluster_mass : array, shape (N_times, N_trials)
       Mass of living stars in non-disrupted clusters at each time in
       each trial
    num_clusters : array, shape (N_times, N_trials), dtype ulonglong
       Number of non-disrupted clusters present at each time in each
       trial
    num_dis_clusters : array, shape (N_times, N_trials), dtype ulonglong
       Number of disrupted clusters present at each time in each trial
    num_fld_stars : array, shape (N_times, N_trials), dtype ulonglong
       Number of living field stars (excluding those in disrupted 
       clusters and those being treated non-stochastically) present at
       each time in each trial

    (Only present if an integrated_spec file is found)
    wl : array
       wavelengths, in Angstrom
    spec : array, shape (N_wavelength, N_times, N_trials)
       specific luminosity at each wavelength and each time for each
       trial, in erg/s/A

    (Only present if an integrated_phot file is found)
    filter_names : list of string
       a list giving the name for each filter
    filter_units : list of string
       a list giving the units for each filter
    filter_wl_cen : list
       central wavelength of each filter; this is set to None for the
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
    """

    # Read properties
    try:
        prop = read_integrated_prop(model_name, output_dir, fmt,
                                    verbose, read_info)
        if read_info is not None:
            read_info['prop_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        prop = None

    # Read spectra
    try:
        spec = read_integrated_spec(model_name, output_dir, fmt,
                                    verbose, read_info)
        if read_info is not None:
            read_info['spec_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        spec = None

    # Read photometry
    try:
        phot = read_integrated_phot(model_name, output_dir, fmt, 
                                    nofilterdata, photsystem, 
                                    verbose, read_info)
        if read_info is not None:
            read_info['phot_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        phot = None

    # Read cloudy spectra
    try:
        cloudyspec = read_integrated_cloudyspec(model_name, output_dir, fmt,
                                                verbose, read_info)
        if read_info is not None:
            read_info['cloudyspec_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        cloudyspec = None


    # Build the output
    out_fields = ['time']
    if prop is not None:
        out_data = [prop.time]
    elif spec is not None:
        out_data = [spec.time]
    elif phot is not None:
        out_data = [phot.time]
    elif cloudyspec is not None:
        out_data = [cloudyspec.time]
    else:
        raise IOError("unable to open any integrated files for run " +
                      model_name)
    if prop is not None:
        out_fields = out_fields + list(prop._fields[1:])
        out_data = out_data + list(prop[1:])
    if spec is not None:
        out_fields = out_fields + list(spec._fields[1:])
        out_data = out_data + list(spec[1:])
    if phot is not None:
        out_fields = out_fields + list(phot._fields[1:])
        out_data = out_data + list(phot[1:])
    if cloudyspec is not None:
        out_fields = out_fields + list(cloudyspec._fields[1:])
        out_data = out_data + list(cloudyspec[1:])
    out_type = namedtuple('integrated_data', out_fields)
    out = out_type._make(out_data)

    # Return data
    return out
                          
