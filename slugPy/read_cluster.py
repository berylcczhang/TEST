"""
Function to read all cluster data for a SLUG2 run.
"""

from collections import namedtuple
from read_cluster_prop import read_cluster_prop
from read_cluster_phot import read_cluster_phot
from read_cluster_spec import read_cluster_spec

def read_cluster(model_name, output_dir=None, asciionly=False,
                 binonly=False, nofilterdata=False, 
                 photsystem=None, verbose=False):
    """
    Function to read all cluster data for a SLUG2 run.

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

    (Always present)
    id : array, dtype uint
       unique ID of cluster
    trial: array, dtype uint
       which trial was this cluster part of
    time : array
       time at which cluster's properties are being evaluated

    (Present if the run being read contains a cluster_prop file)
    form_time : array
       time when cluster formed
    lifetime : array
       time at which cluster will disrupt
    target_mass : array
       target cluster mass
    actual_mass : array
       actual mass at formation
    live_mass : array
       mass of currently living stars
    num_star : array, dtype ulonglong
       number of living stars in cluster being treated stochastically
    max_star_mass : array
       mass of most massive living star in cluster

    (Present if the run being read contains a cluster_spec file)
    wl : array
       wavelength, in Angstrom
    spec : array, shape (N_cluster, N_wavelength)
       specific luminosity of each cluster at each wavelength, in erg/s/A

    (Present if the run being read contains a cluster_phot file)
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
    phot : array, shape (N_cluster, N_filter)
       photometric value in each filter for each cluster; units are as
       indicated in the units field
       
    Raises
    ------
    IOError, if no photometry file can be opened
    ValueError, if photsystem is set to an unknown values
    """

    # Read properties
    try:
        prop = read_cluster_prop(model_name, output_dir, asciionly,
                                 binonly, verbose)
    except IOError:
        prop = None

    # Read spectra
    try:
        spec = read_cluster_spec(model_name, output_dir, asciionly,
                                 binonly, verbose)
    except IOError:
        spec = None

    # Read photometry
    try:
        phot = read_cluster_phot(model_name, output_dir, asciionly,
                                 binonly, nofilterdata, photsystem,
                                 verbose)
    except IOError:
        phot = None

    # Build the output
    out_fields = ['id', 'trial', 'time']
    if prop is not None:
        out_data = [prop.id, prop.trial, prop.time]
    elif spec is not None:
        out_data = [spec.id, spec.trial, spec.time]
    elif phot is not None:
        out_data = [phot.id, phot.trial, phot.time]
    else:
        raise IOError("unable to open any integrated files for run " +
                      model_name)
    if prop is not None:
        out_fields = out_fields + list(prop._fields[3:])
        out_data = out_data + list(prop[3:])
    if spec is not None:
        out_fields = out_fields + list(spec._fields[3:])
        out_data = out_data + list(spec[3:])
    if phot is not None:
        out_fields = out_fields + list(phot._fields[3:])
        out_data = out_data + list(phot[3:])
    out_type = namedtuple('cluster_data', out_fields)
    out = out_type._make(out_data)

    # Return data
    return out

