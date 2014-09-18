"""
Function to combine integrated data from multiple SLUG2 runs, treating
each input run as a separate set of trials.
"""

from collections import namedtuple
import numpy as np

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

    # Safety check: make sure all input data objects have the same
    # fields
    for i in range(len(data)-1):
        if data[i]._fields != data[i+1]._fields:
            raise ValueError("input data must have identical fields")

    # Combine fields
    new_fields = []
    single_fields = ['time', 'target_mass', 'wl', 'filter_names',
                     'filter_units', 'filter_wl', 'filter_wl_eff',
                     'filter_response', 'filter_beta', 'filter_wl_c']
    for i, f in enumerate(data[0]._fields):

        # For the following fields we just need one copy
        if f in single_fields:
            new_fields.append(data[0][i])

        # All other fields are numpy arrays that just get combined
        # along their last axis
        else:
            new_fields.append(data[0][i])
            for j in range(1, len(data)):
                new_fields[-1] \
                    = np.append(new_fields[-1], data[j][i],
                                axis=data[0][i].ndim-1)

    # Create new output object
    out_type = namedtuple('integrated_data', data[0]._fields)
    combined_data = out_type._make(new_fields)

    # Return
    return combined_data
