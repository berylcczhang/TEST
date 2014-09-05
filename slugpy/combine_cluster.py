"""
Function to combine cluster data from multiple SLUG2 runs, treating
each input run as a separate set of trials. Note that trial and unique
ID numbers are not preserved by this operation.
"""

from collections import namedtuple
import numpy as np

def combine_cluster(data):
    """
    Function to combine cluster data from multiple SLUG2 runs,
    treating each input run as a separate set of trials. Trial and
    cluster unique ID numbers are altered as necessary to avoid
    duplication between the merged data sets.

    Parameters
    ----------
    data : list_like
       A list containing the cluster data for each run, as
       returned by read_cluster

    Returns
    -------
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
    for i, f in enumerate(data[0]._fields):

        if f == 'id':
            # ID field requires special handling
            cluster_id = data[0].id
            for j in range(1,len(data)):
                cluster_id = np.append(cluster_id, 
                                       data[j].id+np.amax(cluster_id))
            new_fields.append(cluster_id)

        elif f == 'trial':
            # Trial field requires special handling; note that this is
            # slightly different from the ID field, because the ID
            # field is 1 offset and the trial number is 0 offset
            trial = data[0].trial
            for j in range(1,len(data)):
                trial = np.append(trial,
                                  data[j].trial+np.amax(trial)+1)
            new_fields.append(trial)

        # For the following fields we just need one copy
        elif f == 'wl':
            new_fields.append(data[0].wl)
        elif f == 'filter_names':
            new_fields.append(data[0].filter_names)
        elif f == 'filter_units':
            new_fields.append(data[0].filter_units)
        elif f == 'filter_wl_eff':
            new_fields.append(data[0].filter_wl_eff)
        elif f == 'filter_wl':
            new_fields.append(data[0].filter_wl)
        elif f == 'filter_response':
            new_fields.append(data[0].filter_response)

        # All other fields can just be concatenated
        else:
            joined_field = np.concatenate([d[i] for d in data])
            new_fields.append(joined_field)

    # Create new output object
    out_type = namedtuple('cluster_data', data[0]._fields)
    combined_data = out_type._make(new_fields)

    # Return
    return combined_data
