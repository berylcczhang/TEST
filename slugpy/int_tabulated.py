"""
Code to integrate a tabulated function. The method is a clone of the
IDL int_tabulated function.
"""

import numpy as np
from scipy.interpolate import Akima1DInterpolator as akima

def int_tabulated(x, y, sorted=True):
    """
    Return the approximate integral of y(x), where x is a set of data
    points and y is the function to be integrated evaluated at those
    points.

    Parameters
    ----------
    x : array
       Array of x values
    y : array
       Array of y(x) values
    sorted : bool
       If True, x is assumed to be sorted in increasing order. If it
       is False, it is not assumed to be sorted, and will be sorted
       before integration.

    Returns
    -------
    integral : float
       The integral of y dx
    """

    # Sort
    if not sorted:
        idx = x.argsort()
        x = x[idx]
        y = y[idx]

    # Figure out how many segments we want
    nseg = int(np.ceil((len(x)-1)/4.0)*4)

    # Compute step size and set up interpolation grid; be careful to
    # set the last point exactly so that the interpolator doesn't give
    # us NaN because the last point in the interpolated grid is
    # outside the region of the data due to truncation error
    stepsize = (x[-1] - x[0])/nseg
    xinterp = x[0] + np.arange(nseg+1)*stepsize
    xinterp[-1] = x[-1]

    # Interpolate data onto grid
    finterp = akima(x, y)
    yinterp = finterp(xinterp)

    # Do 5th order Newton-Cotes integration
    total = 2.0 * stepsize * \
            np.sum(7.0*(yinterp[:-1:4] + yinterp[4::4]) +
                   32.0*(yinterp[1::4] + yinterp[3::4]) +
                   12.0*yinterp[2::4]) / 45.0
    return total


def int_tabulated2(x1, f1, x2, f2, sorted=True):
    """
    Return the approximate integral of f1(x) * f2(x), where the
    functions f1 and f2 are tabulated at a set of data points x1 and
    x2.

    Parameters
    ----------
    x1 : array
       Array of x values at which f1 is evaluated
    f1 : array
       Array of f1(x) values
    x2 : array
       Array of x values at which f2 is evaluated
    f2 : array
       Array of f2(x) values
    sorted : bool
       If True, x1 and x2 are assumed to be sorted in increasing
       order. If it is False, they are not assumed to be sorted, and
       will be sorted before integration.

    Returns
    -------
    integral : float
       The integral of f1 f2 dx
    """

    # Sort
    if not sorted:
        idx = x1.argsort()
        x1 = x1[idx]
        f1 = f1[idx]
        idx = x2.argsort()
        x2 = x2[idx]
        f2 = f2[idx]

    # Find the overlapping range
    xmin = max(x1[0], x2[0])
    xmax = min(x1[-1], x2[-1])

    # Construct a new x vector, consisting of all the unique data
    # points in the overlapping range
    x1sub = x1[np.logical_and(x1 >= xmin, x1 <= xmax)]
    x2sub = x2[np.logical_and(x2 >= xmin, x2 <= xmax)]
    x = np.append(x1sub, x2sub)
    x.sort()
    x = np.unique(x)

    # Figure out how many segments we want
    nseg = int(np.ceil((len(x)-1)/4.0)*4)

    # Compute step size and set up interpolation grid
    stepsize = (x[-1] - x[0])/nseg
    xinterp = x[0] + np.arange(nseg+1)*stepsize
    xinterp[-1] = x[-1]    # Avoid problems due to truncation error

    # Interpolate data onto grid
    interp1 = akima(x1, f1)
    f1interp = interp1(xinterp)
    interp2 = akima(x2, f2)
    f2interp = interp2(xinterp)

    # Do 5th order Newton-Cotes integration
    total = 2.0 * stepsize * \
            np.sum(7.0*(f1interp[:-1:4]*f2interp[:-1:4] +
                        f1interp[4::4]*f2interp[4::4]) +
                   32.0*(f1interp[1::4]*f2interp[1::4] + 
                         f1interp[3::4]*f2interp[3::4]) +
                   12.0*f1interp[2::4]*f2interp[2::4]) / 45.0
    return total
