"""
Function to combine cluster data from multiple SLUG2 runs, treating
each input run as a separate set of trials. Note that trial and unique
ID numbers are not preserved by this operation.
"""

from collections import namedtuple
import numpy as np
from warnings import warn

def combine_cluster(data):
    """
    Function to combine cluster data from multiple SLUG2 runs,
    treating each input run as a separate set of trials. Trial and
    cluster unique ID numbers are altered as necessary to avoid
    duplication between the merged data sets.

    Parameters:
       data : list_like
          A list containing the cluster data for each run, as
          returned by read_cluster

    Returns:
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

    # List of fields for which we need only one copy, because they're
    # the same for every cluster
    single_fields = ['wl', 'wl_neb', 'wl_ex', 'wl_neb_ex', 'wl_r',
                     'filter_names', 'filter_units', 'filter_wl',
                     'filter_wl_eff', 'filter_response', 'filter_beta', 
                     'filter_wl_c', 'isotope_name', 'isotope_Z', 'isotope_A',
                     'cloudy_linelabel', 'cloudy_linewl',
                     'cloudy_wl', 'cloudy_filter_names',
                     'cloudy_filter_units', 'cloudy_filter_wl_eff',
                     'cloudy_filter_wl', 'cloudy_filter_response',
                     'cloudy_filter_beta', 'cloudy_filter_wl_c','line_names']

    # Combine fields
    new_fields = []
    for f in s:

        if f == 'id':
            # ID field requires special handling
            cluster_id = data[0].id
            for j in range(1,len(data)):
                if len(cluster_id) > 0:
                    cluster_id = np.append(
                        cluster_id, 
                        data[j].id+np.amax(cluster_id))
                else:
                    cluster_id = np.append(cluster_id, data[j].id)
            new_fields.append(cluster_id)

        elif f == 'trial':
            # Trial field requires special handling; note that this is
            # slightly different from the ID field, because the ID
            # field is 1 offset and the trial number is 0 offset
            trial = data[0].trial
            for j in range(1,len(data)):
                if len(trial) > 0:
                    trial = np.append(
                        trial,
                        data[j].trial+np.amax(trial)+1)
                else:
                    trial = np.append(
                        trial, data[j].trial)
            new_fields.append(trial)

        # For the following fields we just need one copy
        elif f in single_fields:
            new_fields.append(getattr(data[0], f))

        # All other fields can just be concatenated
        else:
            joined_field = np.concatenate([getattr(d, f) for d in data])
            new_fields.append(joined_field)

    # Create new output object
    out_type = namedtuple('cluster_data', s)
    combined_data = out_type._make(new_fields)

    # Return
    return combined_data
