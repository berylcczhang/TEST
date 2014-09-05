"""
Routine to read a cloudy line list output, produced by save last line list
"""

import numpy as np
from collections import namedtuple

def read_cloudy_linelist(filename):
    """
    Reads a cloudy line list output, produced by save last line list

    Parameters
    ----------
    filename : string
       name of the file to be read

    Returns
    -------
    A namedtuple containing the following fields:
    labels : array, dtype 'S4'
       list of line labels
    wl : array
       array of line wavelengths, in Angstrom
    lum : array
       array of line luminosities; this will be in whatever units the
       cloudy output is in
    """

    # Open file
    fp = open(filename, 'r')

    # Burn header lines
    fp.readline()
    fp.readline()

    # Read data
    labels=[]
    wl=[]
    lum=[]
    for line in fp:

        # Skip blank lines
        if len(line.strip()) == 0:
            continue

        # First 4 characters are descriptor
        labels.append(line[:4])

        # Characters 5-10 are wavelength; last character is unit
        wlnum = line[5:10]
        wlunit = line[10]
        if wlunit == 'A':
            multiplier = 1.0
        elif wlunit == 'm':
            multiplier = 1e4
        elif wlunit == ' ':
            multiplier = 1.0
        else:
            raise ValueError("unknown unit in line "+line)
        wl.append(float(wlnum)*multiplier)

        # Quantity after tab is luminosity
        linesplit=line.split('\t')
        lum.append(float(linesplit[1]))

    # Close file
    fp.close()

    # Convert to arrays
    labels = np.array(labels)
    wl = np.array(wl)
    lum = np.array(lum)

    # Return
    out_type = namedtuple('cloudy_linelist', ['label', 'wl', 'lum'])
    return out_type(labels, wl, lum)
