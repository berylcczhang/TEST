"""
Routine to ready cloudy hydrogen conditions and hydrogen ionization
files, produced by 'save last hydrogen conditions' and 'save last
hydrogen ionization'.
"""

import numpy as np
try:
    # Various constants in cgs
    c = physcons.c * 1e2
    a0 = physcons.physical_constants['Bohr radius'][0]*1e2
    alpha = physcons.alpha
    sigma_PI = 2.**9*np.pi**2/(3.0*np.exp(1.0)**4)*alpha*a0**2
except:
    # This exception is here to deal with readthedocs not having scipy
    c = 3.0e10
    sigma_PI = 6.3e-18
fi = 1.1

def read_cloudy_hcon(hcon_file, r0 = 0.0):
    """
    Reads cloudy outputs produce by the 'save last hydrogen
    conditions' and 'save last hydrogen ionization' file, and uses
    these to return various HII region diagnostic parameters.

    Parameters
       hcon_file : str
          Name of hydrogen conditions file to be read
       r0 : float
          Inner radius for the calculation

    Returns
       r1 : float
          outer radius, in cm
       nII : float
          average number density of H nuclei
       Omega : float
          wind parameter, defined as r0^3 / (r1^3 - r0^3)

    Notes
       All averages are averages over the ionized volume, i.e., the
       average is taken with a weighting factor 4 pi r^2 (n_H+/n_H) dV.
    """

    # Load the hydrogen conditions file
    data = np.atleast_2d(np.loadtxt(hcon_file))

    # Form the radii and the weighting factor for volume integrations;
    # be sure to handle pathological case of only 1 or 2 zones of
    # cloudy output
    r = data[:,0] + r0
    if data.shape[0] >= 3:
        dr = 0.5*(data[2:,0] - data[:-2,0])
        dr = np.append(dr, data[-1,0]-data[-2,0])
        dr = np.insert(dr, 0, 0.5*(r[1]+r[0])-r0)
    elif data.shape[0] == 2:
        dr = [ 0.5*(r[1]+r[0])-r0, data[1,0] - data[0,0] ]
    else:
        dr = 1.0   # Doesn't matter in this case, since there's
                   # nothing to weight
    w = r**2*dr*data[:,4]
    wsum = np.sum(w)

    # Get r1 and Omega
    r1 = r[-1]
    Omega = r0**3 / (r1**3 - r0**3)

    # Get average number density
    nII = np.sum(data[:,2]*w)/wsum

    # Return
    return r1, nII, Omega
