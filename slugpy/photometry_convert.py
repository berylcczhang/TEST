"""
Function to convert photometric data between the photometric systems
that SLUG knows about.
"""

import numpy as np
import scipy.constants as physcons

# Constants and units conversions to cm
c = physcons.c*100.0
Angstrom = 1e-8
pc = 3.0856775814671918e18

def photometry_convert(photsystem, phot, units, wl_cen=None,
                       filter_last=False):
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
          different photometric filters
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
                         photsystem)

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

        except NameError:
            raise ValueError("Requested photometric " +
                             "conversion requires " +
                             "filter wavelengths")

        units[i] = target_units

