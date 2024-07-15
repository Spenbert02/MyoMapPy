"""
utility functions for tractography procedures
"""


import numpy as np
from nibabel.streamlines.array_sequence import ArraySequence


def unpackCArray(c_array, lengths, min_num_points=2):
    """from a c-returned list of streamlines (3d array), unpack into and return a nibabel.ArraySequence object containing
    the streamlines

    Parameters
    ----------
    c_array : nd.array[float]
        shape (n, myomappy.tractography.MAX_SL_POINTS, 3) nd.array, where n is the number of streamlines
    lengths : nd.array[int]
        shape (n, ) nd.array, the lengths of each of the n streamline arrays stored in c_array
    min_num_points: int
        the minimum number of points for a streamline. e.g., for min_num_points=1, all length-0
        streamlines in c_array will be ignored
    
    Returns
    -------
    nibabel.ArraySequence
        length-m ArraySequence object, containing the properly truncated streamlines
    """

    ret_arr_seq = ArraySequence()

    for i in range(lengths.shape[0]):
        if lengths[i] < min_num_points:  # current streamline doesn't meet min num points requirement, ignore it
            continue
        curr_np_arr = c_array[i, :lengths[i], :]
        ret_arr_seq.append(curr_np_arr)
    
    return ret_arr_seq
