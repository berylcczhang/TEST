"""
Function to combine integrated data from multiple SLUG2 runs, treating
each input run as a separate set of trials.
"""

from collections import namedtuple
import numpy as np
from warnings import warn

def combine_integrated(data):
    """
    Function to combine integrated data from multiple SLUG2 runs,
    treating each input run as a separate set of trials.

    Parameters
       data : list_like
          A list containing the integrated data for each run, as
          returned by read_integrated

    Returns
       combined_data : namedtuple
          The combined data, in the same format as each object in data
    """

    # Construct list of common fields, and issue warning if fields are
    # not identical
    s = set(data[0]._fields)
    warned = False
    for i in range(len(data)-1):
        if data[i]._fields != data[i+1]._fields and not warned:
            warnstr = "data to be combined does not have identical " + \
                      "fields; only common fields will be kept"
            warn(warnstr)
        s &= set(data[i+1]._fields)

    # Combine fields
    new_fields = []
    single_fields = ['time', 'wl', 'wl_neb', 'wl_ex', 'wl_neb_ex',
                     'filter_names', 'filter_units', 'filter_wl', 
                     'filter_wl_eff', 'filter_response', 'filter_beta', 
                     'filter_wl_c', 'isotope_name', 'isotope_Z', 'isotope_A',
                     'cloudy_linelabel', 'cloudy_linewl',
                     'cloudy_wl', 'cloudy_filter_names',
                     'cloudy_filter_units', 'cloudy_filter_wl_eff',
                     'cloudy_filter_wl', 'cloudy_filter_response',
                     'cloudy_filter_beta', 'cloudy_filter_wl_c']
    for f in s:

        # For the following fields we just need one copy
        if f in single_fields:
            new_fields.append(getattr(data[0], f))

        # All other fields are numpy arrays that just get combined
        # along their last axis
        else:
            new_fields.append(getattr(data[0], f))
            for j in range(1, len(data)):
                new_fields[-1] \
                    = np.append(new_fields[-1], getattr(data[j], f),
                                axis=getattr(data[0], f).ndim-1)

    # Create new output object
    out_type = namedtuple('integrated_data', s)
    combined_data = out_type._make(new_fields)

    # Return
    return combined_data
