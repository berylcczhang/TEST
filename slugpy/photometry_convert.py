"""
Function to convert photometric data between the photometric systems
that SLUG knows about.
"""

import os
import os.path as osp
import numpy as np
import scipy.constants as physcons
from .int_tabulated import int_tabulated
from .int_tabulated import int_tabulated2
from .read_filter import read_filter

# Constants and units conversions to cm
try:
    c = physcons.c*100.0
except:
    # This exception is here to deal with readthedocs not having scipy
    c = 3.0e10
Angstrom = 1e-8
pc = 3.0856775814671918e18

def photometry_convert(photsystem, phot, units, wl_cen=None,
                       filter_last=False, filter_names=None,
                       filter_dir=None):
    """
    Function to convert photometric data between photometric systems.

    Parameters
       photsystem : string
          The photometric system to which to convert. Allowable values
          are 'L_nu', 'L_lambda', 'AB', 'STMAG', and 'Vega',
          corresponding to the options defined in the SLUG code. If this
          is set and the conversion requested involves a conversion from
          a wavelength-based system to a frequency-based one, wl_cen must
          not be None.
       phot : array
          array of photometric data; if the array has more than one
          dimension, the first dimension is assumed to represent the
          different photometric filters (unless filter_last is True,
          in which case the last dimension is represents the array of
          filters)
       units : iterable of strings
          iterable listing the units of the input photometric data. On
          return, strings will be changed to the units of the new system.
       wl_cen : array
          central wavelengths of the filters, in Angstrom; can be left as
          None if the requested conversion doesn't require going between
          wavelength- and frequency-based systems.
       filter_last : bool
          If the input data have more than one dimension, by default it
          is assumed that the first dimension contains values for the
          different photometric filters. If this keyword is set to True,
          it will instead be assumed that the last dimension contains the
          values for the different filters.
       filter_names : iterable of strings
          Names of all filters, used to read the filter response
          functions from disk; only needed for conversions to and from
          Vega magnitudes, and ignored otherwise
       filter_dir : string
          Directory where the filter data files can be found. If left as
          None, filters will be looked for in the $SLUG_DIR/lib/filters
          directory. This parameter is used only for conversions to
          and from Vega magnitudes.

    Returns
       Nothing

    Raises
       ValueError, if wl_cen is None but the requested conversion
       requires going between wavelength- and frequency-based systems
    """

    # Set target units based on requested conversion
    if photsystem == 'L_lambda':
        target_units = 'erg/s/A'
    elif photsystem == 'L_nu':
        target_units = 'erg/s/Hz'
    elif photsystem == 'AB':
        target_units = 'AB mag'
    elif photsystem == 'STMAG':
        target_units = 'ST mag'
    elif photsystem == 'Vega':
        target_units = 'Vega mag'
    else:
        raise ValueError('Unknown photometric system ' +
                         repr(photsystem))

    # If Vega magnitudes is either the current unit or the desination
    # unit, compute the magnitude of Vega in each of our filters
    if photsystem == 'Vega' or 'Vega mag' in units:

        # Make sure we have filter names
        if filter_names is None:
            raise ValueError("must set filter_names keyword for "
                             "conversions to or from Vega mag")

        # Look in SLUG_DIR for Vega spectrum
        from astropy.io import ascii
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
                raise IOError(1, "could not find Vega spectrum "+
                              "file A0V_KURUCZ_92.SED in "+
                              "SLUG_DIR/lib/atmospheres or "+
                              "in current directory")

        # Record Vega data
        vega_logwl = np.log(np.array(vegadata['col1']))
        vega_spec = np.array(vegadata['col2'])
        vega_F_nu = Angstrom * vega_spec * np.exp(2.0*vega_logwl) / c

        # Get the convolution of Vega's spectrum with the Johnson V
        # filter, and convert that to magnitudes
        fdat = read_filter('Johnson_V', filter_dir=filter_dir)
        fnorm = int_tabulated(np.log(fdat.wl), fdat.response)
        vega_Fbar_nu \
            = int_tabulated2(vega_logwl, vega_F_nu,
                             np.log(fdat.wl), fdat.response) / fnorm
        vega_mag_V = -2.5*np.log10(vega_Fbar_nu)

        # Now loop over other filters, and compute the magnitude of
        # Vega in each one, normalizing to give a magnitude of 0.02 in
        # V
        vega_mag = np.zeros(len(filter_names))
        for i in range(len(filter_names)):
            # Skip special values
            if (units[i] == 'Lsun') or (units[i] == 'phot/s'):
                continue
            fdat = read_filter(filter_names[i], filter_dir=filter_dir)
            fnorm = int_tabulated(np.log(fdat.wl), fdat.response)
            vega_Fbar_nu_filter \
                = int_tabulated2(vega_logwl, vega_F_nu,
                                 np.log(fdat.wl), 
                                 fdat.response) / fnorm
            vega_mag_filter = -2.5*np.log10(vega_Fbar_nu_filter)
            vega_mag[i] = vega_mag_filter - vega_mag_V + 0.02

    # If our input data are in Vega mag, convert them to AB mag first
    for i in range(len(units)):
        if units[i] == 'Vega mag':
            if filter_last:
                phot[:,i] = phot[:,i] + vega_mag[i]
            else:
                phot[i,:] = phot[i,:] + vega_mag[i]
            units[i] = 'AB mag'

    # Loop over photometric values
    for i in range(len(units)):

        # If units are Lsun or phot/s, this is a special value
        # that should not be converted. Also do nothing if target
        # units match current units.
        if (units[i] == 'Lsun') or (units[i] == 'phot/s') or \
           (units[i] == target_units):
            continue

        # Do unit conversion
        try:

            # L_lambda to other units
            if (units[i] == 'erg/s/A') and \
               (target_units == 'erg/s/Hz'):
                if filter_last:
                    phot[:,i] = phot[:,i] * wl_cen[i]**2 * Angstrom / c
                else:
                    phot[i,:] = phot[i,:] * wl_cen[i]**2 * Angstrom / c
            elif (units[i] == 'erg/s/A') and \
                 (target_units == 'AB mag'):
                if filter_last:
                    L_nu = phot[:,i] * wl_cen[i]**2 * Angstrom / c
                    F_nu = L_nu / (4.0*np.pi*(10.0*pc)**2)
                    phot[:,i] = -2.5*np.log10(F_nu) - 48.6
                else:
                    L_nu = phot[i,:] * wl_cen[i]**2 * Angstrom / c
                    F_nu = L_nu / (4.0*np.pi*(10.0*pc)**2)
                    phot[i,:] = -2.5*np.log10(F_nu) - 48.6
            elif (units[i] == 'erg/s/A') and \
                 (target_units == 'Vega mag'):
                if filter_last:
                    L_nu = phot[:,i] * wl_cen[i]**2 * Angstrom / c
                    F_nu = L_nu / (4.0*np.pi*(10.0*pc)**2)
                    phot[:,i] = -2.5*np.log10(F_nu) - 48.6 - vega_mag[i]
                else:
                    L_nu = phot[i,:] * wl_cen[i]**2 * Angstrom / c
                    F_nu = L_nu / (4.0*np.pi*(10.0*pc)**2)
                    phot[i,:] = -2.5*np.log10(F_nu) - 48.6 - vega_mag[i]
            elif (units[i] == 'erg/s/A') and \
                 (target_units == 'ST mag'):
                if filter_last:
                    F_lambda = phot[:,i] / (4.0*np.pi*(10.0*pc)**2)
                    phot[:,i] = -2.5*np.log10(F_lambda) - 21.1
                else:
                    F_lambda = phot[i,:] / (4.0*np.pi*(10.0*pc)**2)
                    phot[i,:] = -2.5*np.log10(F_lambda) - 21.1

            # L_nu to other units
            elif (units[i] == 'erg/s/Hz') and \
                 (target_units == 'erg/s/A'):
                if filter_last:
                    phot[:,i] = phot[:,i] * c / \
                                (wl_cen[i]**2 * Angstrom)
                else:
                    phot[i,:] = phot[i,:] * c / \
                                (wl_cen[i]**2 * Angstrom)
            elif (units[i] == 'erg/s/Hz') and \
                 (target_units == 'AB mag'):
                if filter_last:
                    F_nu = phot[:,i] / (4.0*np.pi*(10.0*pc)**2)
                    phot[:,i] = -2.5*np.log10(F_nu) - 48.6
                else:
                    F_nu = phot[i,:] / (4.0*np.pi*(10.0*pc)**2)
                    phot[i,:] = -2.5*np.log10(F_nu) - 48.6
            elif (units[i] == 'erg/s/Hz') and \
                 (target_units == 'Vega mag'):
                if filter_last:
                    F_nu = phot[:,i] / (4.0*np.pi*(10.0*pc)**2)
                    phot[:,i] = -2.5*np.log10(F_nu) - 48.6 - vega_mag[i]
                else:
                    F_nu = phot[i,:] / (4.0*np.pi*(10.0*pc)**2)
                    phot[i,:] = -2.5*np.log10(F_nu) - 48.6 - vega_mag[i]
            elif (units[i] == 'erg/s/Hz') and \
                 (target_units == 'ST mag'):
                if filter_last:
                    L_lambda = phot[:,i] * c / (wl_cen[i]**2 * Angstrom)
                    F_lambda = L_lambda / (4.0*np.pi*(10.0*pc)**2)
                    phot[:,i] = -2.5*np.log10(F_lambda) - 21.1
                else:
                    L_lambda = phot[i,:] * c / (wl_cen[i]**2 * Angstrom)
                    F_lambda = L_lambda / (4.0*np.pi*(10.0*pc)**2)
                    phot[i,:] = -2.5*np.log10(F_lambda) - 21.1

            # AB mag to other units
            elif (units[i] == 'AB mag') and \
                 (target_units == 'erg/s/A'):
                if filter_last:
                    F_nu = 10.0**(-(phot[:,i] + 48.6)/2.5)
                    F_lambda = F_nu * c / (wl_cen[i]**2 * Angstrom)
                    phot[:,i] = F_lambda * 4.0*np.pi*(10.0*pc)**2
                else:
                    F_nu = 10.0**(-(phot[i,:] + 48.6)/2.5)
                    F_lambda = F_nu * c / (wl_cen[i]**2 * Angstrom)
                    phot[i,:] = F_lambda * 4.0*np.pi*(10.0*pc)**2
            elif (units[i] == 'AB mag') and \
                 (target_units == 'erg/s/Hz'):
                if filter_last:
                    F_nu = 10.0**(-(phot[:,i] + 48.6)/2.5)
                    phot[:,i] = F_nu * 4.0*np.pi*(10.0*pc)**2
                else:
                    F_nu = 10.0**(-(phot[i,:] + 48.6)/2.5)
                    phot[i,:] = F_nu * 4.0*np.pi*(10.0*pc)**2
            elif (units[i] == 'AB mag') and \
                 (target_units == 'ST mag'):
                if filter_last:
                    F_nu = 10.0**(-(phot[:,i] + 48.6)/2.5)
                    F_lambda = F_nu * c / (wl_cen[i]**2 * Angstrom)
                    phot[:,i] = -2.5*np.log10(F_lambda) - 21.1
                else:
                    F_nu = 10.0**(-(phot[i,:] + 48.6)/2.5)
                    F_lambda = F_nu * c / (wl_cen[i]**2 * Angstrom)
                    phot[i,:] = -2.5*np.log10(F_lambda) - 21.1
            elif (units[i] == 'AB mag') and \
                 (target_units == 'Vega mag'):
                if filter_last:
                    phot[:,i] = phot[:,i] - vega_mag[i]
                else:
                    phot[i,:] = phot[i,:] - vega_mag[i]

            # ST mag to other units
            elif (units[i] == 'ST mag') and \
                 (target_units == 'erg/s/A'):
                if filter_last:
                    F_lambda = 10.0**(-(phot[:,i] + 21.1)/2.5)
                    phot[:,i] = F_lambda * 4.0*np.pi*(10.0*pc)**2
                else:
                    F_lambda = 10.0**(-(phot[i,:] + 21.1)/2.5)
                    phot[i,:] = F_lambda * 4.0*np.pi*(10.0*pc)**2
            elif (units[i] == 'ST mag') and \
                 (target_units == 'erg/s/Hz'):
                if filter_last:
                    F_lambda = 10.0**(-(phot[:,i] + 21.1)/2.5)
                    F_nu = F_lambda * wl_cen[i]**2 * Angstrom / c
                    phot[:,i] = F_nu * 4.0*np.pi*(10.0*pc)**2
                else:
                    F_lambda = 10.0**(-(phot[i,:] + 21.1)/2.5)
                    F_nu = F_lambda * wl_cen[i]**2 * Angstrom / c
                    phot[i,:] = F_nu * 4.0*np.pi*(10.0*pc)**2
            elif (units[i] == 'ST mag') and \
                 (target_units == 'AB mag'):
                if filter_last:
                    F_lambda = 10.0**(-(phot[:,i] + 21.1)/2.5)
                    F_nu = F_lambda * wl_cen[i]**2 * Angstrom / c
                    phot[:,i] = -2.5*np.log10(F_nu) - 48.6
                else:
                    F_lambda = 10.0**(-(phot[i,:] + 21.1)/2.5)
                    F_nu = F_lambda * wl_cen[i]**2 * Angstrom / c
                    phot[i,:] = -2.5*np.log10(F_nu) - 48.6
            elif (units[i] == 'ST mag') and \
                 (target_units == 'Vega mag'):
                if filter_last:
                    F_lambda = 10.0**(-(phot[:,i] + 21.1)/2.5)
                    F_nu = F_lambda * wl_cen[i]**2 * Angstrom / c
                    phot[:,i] = -2.5*np.log10(F_nu) - 48.6 - vega_mag[i]
                else:
                    F_lambda = 10.0**(-(phot[i,:] + 21.1)/2.5)
                    F_nu = F_lambda * wl_cen[i]**2 * Angstrom / c
                    phot[i,:] = -2.5*np.log10(F_nu) - 48.6 - vega_mag[i]

        except NameError:
            raise ValueError("Requested photometric " +
                             "conversion requires " +
                             "filter wavelengths")

        # Store new units
        units[i] = target_units

