"""
Routine to read a cloudy continuum output, produced by save last continuum
"""

import numpy as np
from collections import namedtuple
import scipy.constants as physcons
try:
    c = physcons.c * 1e2   # c in CGS units
except:
    # This exception is here to deal with readthedocs not having scipy
    c = 3.0e10

def read_cloudy_continuum(filename, r0=None):
    """
    Reads a cloudy continuum output, produced by save last continuum

    Parameters
       filename : string
          name of the file to be read
       r0 : float
          inner radius, in cm; if included, the quantities returned will
          be total energies instead of energy emission rates instead of
          rates per unit area

    Returns
       A namedtuple containing the following fields:

       wl : array
          wavelengths in Angstrom
       incident : array
          incident radiation field intensity
    """

    # Open file
    fp = open(filename, 'r')

    # Burn the header line
    fp.readline()

    # Prepare storage
    en = []
    emission = []
    line_label = []
    cont_label = []
    line_density = []

    # Read the data
    for line in fp:
        # Skip comment lines that start with #
        if line.strip()[0] == '#':
            continue
        linesplit = line.split()
        en.append(float(linesplit[0]))
        incident = float(linesplit[1])
        transmitted = float(linesplit[2])
        emitted = float(linesplit[3])
        transmitted_plus_emitted = float(linesplit[4])
        emission.append([incident, transmitted, emitted,
                         transmitted_plus_emitted])
    fp.close()

    # Convert to arrays
    en = np.array(en)
    emission = np.transpose(np.array(emission))

    # Convert units to match convension in SLUG -- emission is
    # L_lambda, units of erg / s / Ang, flux is F_lambda, units of erg
    # / cm^2 / s / Ang; cloudy gives us nu L_nu / 4 pi r_0^2, units of
    # erg / s / cm^2
    nu = en * physcons.Rydberg*physcons.c  # Rydbergs to Hz
    F_nu = emission / nu                   # Convert nu L_nu to L_nu
    F_lambda = F_nu * nu**2/c / 1e8        # Convert L_nu to L_lambda,
                                           # and units to erg / s / Ang
    wl = 1e8 * c/nu[::-1]                  # Wavelength in Angstrom
    F_lambda = np.copy(F_lambda[:,::-1])   # Reverse order
    if r0 is not None:
        L_lambda = F_lambda * 4.0*np.pi*r0**2
        out_type = namedtuple('cloudy_continuum',
                              ['wl', 'L_lambda'])
        out = out_type(wl, L_lambda)
    else:
        out_type = namedtuple('cloudy_continuum',
                              ['wl', 'F_lambda'])
        out = out_type(wl, F_lambda)

    # Return output
    return out
