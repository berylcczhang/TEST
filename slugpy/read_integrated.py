"""
Function to read all integrated data for a SLUG2 run.
"""

from collections import namedtuple
from .read_integrated_prop import read_integrated_prop
from .read_integrated_phot import read_integrated_phot
from .read_integrated_spec import read_integrated_spec
from .read_integrated_yield import read_integrated_yield
from .cloudy.read_integrated_cloudyparams import read_integrated_cloudyparams
from .cloudy.read_integrated_cloudyphot import read_integrated_cloudyphot
from .cloudy.read_integrated_cloudylines import read_integrated_cloudylines
from .cloudy.read_integrated_cloudyspec import read_integrated_cloudyspec

def read_integrated(model_name, output_dir=None, fmt=None,
                    nofilterdata=False, photsystem=None, 
                    verbose=False, read_info=None,
                    no_stellar_mass=False):
    """
    Function to read all integrated light data for a SLUG2 run.

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
          On return, this dict will contain the keys 'prop_name',
          'phot_name', 'spec_name', 'cloudyspec_name', 'cloudylines_name'
          and 'format', giving the names of the files read and the format
          they were in; 'format' will be one of 'ascii', 'binary', or
          'fits'. If one of the files is not present, the corresponding
          _name key will be omitted from the dict.
       no_stellar_mass : bool
          Prior to 7/15, output files did not contain the stellar_mass
          field; this can be detected automatically for ASCII and FITS
          formats, but not for binary format; if True, this specifies
          that the binary file being read does not contain a
          stellar_mass field; it has no effect for ASCII or FITS files

    Returns
       A namedtuple containing the following fields:

       (Always present)

       time: array
          Times at which data are output

       (Only present if an integrated_prop file is found)

       target_mass : array, shape (N_times, N_trials)
          Target stellar mass at each time in each trial
       actual_mass : array, shape (N_times, N_trials)
          Actual mass of stars created up to each time in each trial
       live_mass : array, shape (N_times, N_trials)
          Mass of currently-alive stars at each time in each trial
       stellar_mass : array
          mass of all stars, living and stellar remnants
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
       wl_neb : array
          wavelength for the nebular spectrum, in Angstrom (present
          only if SLUG was run with nebular emission enabled)
       spec_neb : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity at each wavelength and each time for each
          trial, including emission and absorption by the HII region,
          in erg/s/A (present only if SLUG was run with nebular
          emission enabled)
       wl_ex : array
          wavelength for the extincted spectrum, in Angstrom (present
          only if SLUG was run with extinction enabled)
       spec_ex : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity at each wavelength in wl_ex and each
          time for each trial after extinction has been applied, in
          erg/s/A (present only if SLUG was run with extinction
          enabled)
       wl_neb_ex : array
          wavelength for the extincted spectrum with nebular emission,
          in Angstrom (present only if SLUG was run with both nebular
          emission and extinction enabled)
       spec_neb_ex : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity at each wavelength in wl_ex and each
          time for each trial including emission and absorption by the
          HII region, after extinction has been applied, in erg/s/A
          (present only if SLUG was run with nebular emission and
          extinction both enabled)

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

       If extinction is enabled, phot_ex will contain photometry  
       after extinction has been applied.     

       (Only present if an integrate_yield file is found)

       isotope_name : array of strings
          Atomic symbols of the isotopes included in the yield table
       isotope_Z : array of int
          Atomic numbers of the isotopes included in the yield table
       isotope_A : array of int
          Atomic mass number of the isotopes included in the yield table
       yld : array
          Yield of each isotope, defined as the instantaneous amount
          produced up to that time; for unstable isotopes, this
          includes the effects of decay since production

       (Only present if an integrated_cloudyspec file is found)

       cloudy_wl : array
          wavelength, in Angstrom
       cloudy_inc : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity of the stellar radiation field at each
          wavelength and each time for each trial, in erg/s/A
       cloudy_trans : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity of the stellar radiation field after it has
          passed through the HII region, at each wavelength and each time
          for each trial, in erg/s/A
       cloudy_emit : array, shape (N_wavelength, N_times, N_trials)
          specific luminosity of the radiation field emitted by the HII
          region, at each wavelength and each time for each trial, in
          erg/s/A
       cloudy_trans_emit : array, shape (N_wavelength, N_times, N_trials)
          the sum of emitted and transmitted; this is what would be seen
          by an observer looking at both the star cluster and its nebula

       (Only present if an integrated_cloudylines file is found)

       cloudy_linelabel : array, dtype='S4', shape (N_lines)
          labels for the lines, following cloudy's 4 character line label
          notation
       cloudy_linewl : array, shape (N_lines)
          rest wavelength for each line, in Angstrom
       cloudy_linelum : array, shape (N_lines, N_times, N_trials)
          luminosity of each line at each time for each trial, in erg/s

       (Only present if an integrated_cloudyphot file is found)

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

       (Only present if an integrated_cloudyparams file is found)

       cloudy_hden : array
          number density of H nuclei at the inner edge of the ionized
          region simulated by cloudy
       cloudy_r0 : array
          inner radius of the ionized region simulated by cloudy
       cloudy_rS : array
          outer radius of the ionized region simulated by cloudy (approximate!)
       cloudy_QH0 : array
          ionizing luminosity used in the cloudy computation
       cloudy_covFac : array
          covering factor assumed in the cloudy computation; only a
          fraction covFac of the ionizing photons are assumed to
          produce emission within the HII region, while the remainder
          are assumed to escape
       cloudy_U : array
          volume-averaged ionization parameter of the HII region
          simulated by cloudy; note that this value is approximate,
          not exact, and the approximation can be very poor if
          radiation pressure effects are significant
       cloudy_Omega : array
          Yeh & Matzner (2012) wind parameter for the HII region
          simulated by cloudy; as with U, this value is approximate,
          and the approximation is valid only if radiation pressure
          effects are small
    """

    # Read properties
    try:
        prop = read_integrated_prop(model_name, output_dir, fmt=fmt,
                                    verbose=verbose, 
                                    read_info=read_info,
                                    no_stellar_mass=no_stellar_mass)
        if read_info is not None:
            read_info['prop_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        prop = None

    # Read spectra
    try:
        spec = read_integrated_spec(model_name, output_dir, fmt=fmt,
                                    verbose=verbose, 
                                    read_info=read_info)
        if read_info is not None:
            read_info['spec_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        spec = None

    # Read photometry
    try:
        phot = read_integrated_phot(model_name, output_dir, fmt=fmt, 
                                    nofilterdata=nofilterdata,
                                    photsystem=photsystem, 
                                    verbose=verbose,
                                    read_info=read_info)
        if read_info is not None:
            read_info['phot_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        phot = None

    # Read cloudy spectra
    try:
        cloudyspec = read_integrated_cloudyspec(model_name, output_dir, 
                                                fmt=fmt,
                                                verbose=verbose, 
                                                read_info=read_info)
        if read_info is not None:
            read_info['cloudyspec_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        cloudyspec = None

    # Read cloudy lines
    try:
        cloudylines \
            = read_integrated_cloudylines(model_name, output_dir,
                                          fmt=fmt,
                                          verbose=verbose,
                                          read_info=read_info)
        if read_info is not None:
            read_info['cloudylines_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        cloudylines = None

    # Read cloudy photometry
    try:
        cloudyphot \
            = read_integrated_cloudyphot(model_name, output_dir, 
                                         fmt=fmt,
                                         nofilterdata=nofilterdata,
                                         photsystem=photsystem, 
                                         verbose=verbose,
                                         read_info=read_info)
        if read_info is not None:
            read_info['cloudyphot_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        cloudyphot = None

    # Read cloudy parameters
    try:
        cloudyparams \
            = read_integrated_cloudyparams(
                model_name, output_dir, fmt=fmt,
                verbose=verbose, read_info=read_info)
        if read_info is not None:
            read_info['cloudylines_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        cloudyparams = None

    # Read yield data
    try:
        yld = read_integrated_yield(model_name, output_dir=output_dir,
                                    fmt=fmt, verbose=verbose,
                                    read_info=read_info)
        if read_info is not None:
            read_info['yield_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        yld = None

    # Build the output
    out_fields = ['time']
    if prop is not None:
        out_data = [prop.time]
    elif spec is not None:
        out_data = [spec.time]
    elif phot is not None:
        out_data = [phot.time]
    elif yld is not None:
        out_data = [yld.time]
    elif cloudyspec is not None:
        out_data = [cloudyspec.time]
    elif cloudylines is not None:
        out_data = [cloudylines.time]
    elif cloudyphot is not None:
        out_data = [cloudyphot.time]
    elif cloudyparams is not None:
        out_data = [cloudyparams.time]
    else:
        raise IOError(1,
                      "unable to open any integrated files for run " +
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
    if yld is not None:
        out_fields = out_fields + list(yld._fields[1:])
        out_data = out_data + list(yld[1:])
    if cloudyspec is not None:
        out_fields = out_fields + list(cloudyspec._fields[1:])
        out_data = out_data + list(cloudyspec[1:])
    if cloudylines is not None:
        out_fields = out_fields + list(cloudylines._fields[1:])
        out_data = out_data + list(cloudylines[1:])
    if cloudyphot is not None:
        out_fields = out_fields + list(cloudyphot._fields[1:])
        out_data = out_data + list(cloudyphot[1:])
    if cloudyparams is not None:
        out_fields = out_fields + list(cloudyparams._fields[1:])
        out_data = out_data + list(cloudyparams[1:])
    out_type = namedtuple('integrated_data', out_fields)
    out = out_type._make(out_data)

    # Return data
    return out
                          
