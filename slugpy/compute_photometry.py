"""
This function takes an input spectrum and a set of response functions
for photometric filters, and returns the photometry through those
filters.
"""

import os
import os.path as osp
import numpy as np
from .int_tabulated import int_tabulated
from .int_tabulated import int_tabulated2
from .read_filter import read_filter
import scipy.constants as physcons
import errno

# Units and constants
try:
    c = physcons.c*1e2
    h = physcons.h*1e7
except:
    # This exception is here to deal with readthedocs not having scipy
    c = 3.0e10
    h = 6.63e-27
ang = 1e-8
pc = 3.0856775814671918e18
Lsun = 3.846e33

# Threshold wavelengths of H0, He0, He1 ionizing photons
lambda_thresh = { 
    'QH0' : 911.2670505509151,
    'QHe0' : 503.9888,
    'QHe1' : 227.2979
}

def compute_photometry(wl, spec, filtername, photsystem='L_nu', 
                       filter_wl=None, filter_response=None,
                       filter_beta=None, filter_wl_c=None,
                       filter_dir=None):
    """
    This function takes an input spectrum and a set of response
    functions for photometric filters, and returns the photometry
    through those filters.
    
    Parameters
       wl : array
          Wavelength of input spectrum in Angstrom
       spec : array
          Specific luminosity per unit wavelength for input spectrum, in
          erg/s/A
       filtername : string or iterable of strings
          Name or list of names of the filters to be used. Filter names
          can also include the special filters Lbol, QH0, QHe0, and QHe1;
          the values returned for these will be the bolometric luminosity
          (in erg/s) and the photon luminosities (in photons/s) in the H,
          He, and He+ ionizing-continua, respectively.
       photsystem : string
          The photometric system to use for the output. Allowable values
          are 'L_nu', 'L_lambda', 'AB', 'STMAG', and 'Vega',
          corresponding to the options defined in the SLUG code.
       filter_wl : array or iterable of arrays
          Array giving the wavelengths in Angstrom at which the filter is
          response function is given. If this object is an iterable of
          arrays rather than a single array, it is assumed to represent
          the wavelengths for a set of filters. If this is set,
          no data is read from disk. Default behavior is to read the
          filter information from disk.
       filter_response : array or iterable of arrays
          Array giving the filter response function at each wavelenght
          and for each filter in filter_wl. Must be set if filter_wl is
          set, ignored otherwise.
       filter_beta : iterable
          Array-like object containing the index beta for each
          filter. Must be set if filter_wl is set, ignored otherwise.
       filter_wl_c : iterable
          Array-like object containing the pivot wavelength for each
          filter. Must be set if filter_wl is set, ignored otherwise.
       filter_dir : string
          Directory where the filter data files can be found. If left as
          None, filters will be looked for in the $SLUG_DIR/lib/filters
          directory. This parameter is used only if filtername is not
          None.

    Returns
       phot : array
          Photometric values in the requested filters. Units depend on
          the choice of photometric system:
          L_nu --> erg/s/Hz;
          L_lambda --> erg/s/A;
          AB --> absolute AB magnitude;
          STMAG --> absolute ST magnitude;
          Vega --> absolute Vega magnitude;
    """

    # Read filter data if needed
    if filter_wl is None:
        filter_wl = []
        filter_response = []
        filter_beta = []
        filter_wl_c = []
        if not hasattr(filtername, '__iter__'):
            filtername = [filtername]
        for f in filtername:
            filterdata = read_filter(f, filter_dir=filter_dir)
            filter_wl.append(filterdata[1])
            filter_response.append(filterdata[2])
            filter_beta.append(filterdata[3])
            filter_wl_c.append(filterdata[4])

    # Compute normalization for each filter
    filter_norm = []
    filter_logwl = []
    for i in range(len(filter_wl)):
        if filter_wl[i] is None:
            filter_norm.append(None)
            filter_logwl.append(None)
        else:
            filter_logwl.append(np.log(filter_wl[i]))
            if filter_wl_c[i] is None:
                integrand = filter_response[i]
            else:
                integrand = filter_response[i] * \
                            (filter_wl[i]/filter_wl_c[i])**(-filter_beta[i])
            filter_norm.append(int_tabulated(filter_logwl[-1], integrand))

    # Log of wavelength
    logwl = np.log(wl)

    # Loop over filters
    phot = np.zeros(len(filter_wl))
    for i in range(len(filter_wl)):

        # See if this is a special filter
        if filtername is None:
            fname = ''
        else:
            fname = filtername[i]

        if fname == 'Lbol':

            # Bolometric filter, so just integrate and convert to Lsun
            phot[i] = int_tabulated(logwl, spec*wl) / Lsun

        elif fname in lambda_thresh.keys():

            # Ionizing photon luminosity

            # Get log wavelength in the range we want
            logwl_sub = logwl[wl < lambda_thresh[fname]]
            if len(logwl_sub) == 0:
                # No ionizing part to the spectrum
                phot[i] = 0.0
                continue

            # Get integrand, lambda^2 L_lambda / hc
            integrand = spec[wl < lambda_thresh[fname]] * \
                        wl[wl < lambda_thresh[fname]]**2 / (h*c)

            # Integrate
            phot[i] = ang * int_tabulated(logwl_sub, integrand)

        else:

            # All other filters

            # Integrate depending on photometric system
            if photsystem == 'L_nu':

                # L_nu mode
                L_nu = ang * spec * wl**2 / c
                phot[i] = int_tabulated2(logwl, L_nu, 
                                         filter_logwl[i], 
                                         filter_response[i]) / filter_norm[i]

            elif photsystem == 'L_lambda':

                # L_lambda mode
                phot[i] = int_tabulated2(logwl, spec, 
                                         filter_logwl[i], 
                                         filter_response[i]) / filter_norm[i]

            elif photsystem == 'AB' or photsystem == 'Vega':

                # AB mag mode or Vega mag mode; for now treat this as
                # AB in either case, and convert to Vega below if
                # needed
                L_nu = ang * spec * wl**2 / c
                Lbar_nu = int_tabulated2(logwl, L_nu, 
                                         filter_logwl[i], 
                                         filter_response[i]) / filter_norm[i]
                F_nu = Lbar_nu / (4.0*np.pi*(10.0*pc)**2)
                phot[i] = -2.5*np.log10(F_nu) - 48.6

            elif photsystem == 'STMAG':

                # ST mag mode
                Lbar_lambda \
                    = int_tabulated2(logwl, spec, 
                                     filter_logwl[i], 
                                     filter_response[i]) / filter_norm[i]
                F_lambda = Lbar_lambda / (4.0*np.pi*(10.0*pc)**2)
                phot[i] = -2.5*np.log10(F_lambda) - 21.1

            else:

                # Unknown photometry mode, throw error
                raise ValueError("unknown photometric system " +
                                 str(photsystem))


    # For Vega magnitudes, convert from AB to that
    if photsystem == 'Vega':

        from astropy.io import ascii

        # Look in SLUG_DIR for Vega spectrum
        if 'SLUG_DIR' in os.environ:
            slugdir = os.environ['SLUG_DIR']
            fname = osp.join(slugdir, 'lib', 'atmospheres',
                             'A0V_KURUCZ_92.SED')
            try:
                vegadata = ascii.read(fname)
            except IOError:
                vegadata = None

        # Look in cwd if that fails
        if vegadata is None:
            try:
                vegadata = ascii.read('A0V_KURUCZ_92.SED', 'r')
            except IOError:
                raise IOError(errno.ENOENT,
                              "could not find Vega spectrum "+
                              "file A0V_KURUCZ_92.SED in "+
                              "SLUG_DIR/lib/atmospheres or "+
                              "in current directory")

        # Record Vega data
        vega_logwl = np.log(np.array(vegadata['col1']))
        vega_spec = np.array(vegadata['col2'])
        vega_F_nu = ang * vega_spec * np.exp(2.0*vega_logwl) / c

        # Get the convolution of Vega's spectrum with the Johnson V
        # filter, and convert that to magnitudes
        fdat = read_filter('Johnson_V', filter_dir=filter_dir)
        v_norm = int_tabulated(np.log(fdat.wl), fdat.response)
        vega_Fbar_nu \
            = int_tabulated2(vega_logwl, vega_F_nu,
                             np.log(fdat.wl), fdat.response) / v_norm
        vega_mag_offset = -2.5*np.log10(vega_Fbar_nu)

        # Now figure out the magnitude of Vega in all other filters,
        # and use that to convert AB magnitudes to Vega magnitudes
        for i in range(len(filter_wl)):

            # Skip special filters
            if filtername is None:
                fname = ''
            else:
                fname = filtername[i]
            if fname == 'Lbol' or fname in lambda_thresh.keys():
                continue

            # Get magnitudes in this filter before applying offset
            vega_Fbar_nu_filter \
                = int_tabulated2(vega_logwl, vega_F_nu,
                                 filter_logwl[i],
                                 filter_response[i]) / filter_norm[i]
            vega_mag_filter = -2.5*np.log10(vega_Fbar_nu_filter)

            # Compute magnitude of Vega in this band from condition
            # that Vega has magnitude 0.02 in Johnson_V
            vega_mag = vega_mag_filter - vega_mag_offset + 0.02

            # Now convert AB to Vega mag
            phot[i] = phot[i] - vega_mag

    # Return
    return phot
