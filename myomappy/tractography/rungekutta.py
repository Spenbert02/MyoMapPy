"""
Runge Kutta Tractography functions
"""


import ctypes
import numpy as np
import os.path
import platform
from myomappy.tractography import TractographyMethods, MAX_SL_POINTS
from myomappy.tractography._utils import unpackCArray
from numpy.ctypeslib import ndpointer


def rk4Track(dirs, mask, affine, step_size, max_length=None, max_prop_angle=None, seed_density=1):
    """from a direction field and mask, use 4th order Runge-Kutta tractography to
    create and return a set of tracts

    Parameters
    ----------
    dirs : ndarray
        (M, N, P, 3) shape nd.array, the per-voxel direction field to track
    mask : ndarray
        (M, N, P) shape nd.array, the binary voxel mask to track within
    affine : ndarray
        (4, 4) shape nd.array, the voxel-to-world homogenous transformation matrix
    step_size : float, int
        the step size to use during tracking
    max_length : float, int, None
        the (inclusive) maximum tract length. Upon reaching this length, tracts are terminated. If
        None, there is no limit on tract length
    max_prop_angle : float, int, None
        the maximum tract propagation angle (degrees). If this angle is exceeded during tracking, the
        tract is terminated. If None, propagation angle is not a termination criteria.
    seed_density: int
        the seed density to use for tractography. e.g., a seed density of 2 uses 2x2x2 seed points
        per voxel.

    Returns
    -------
    nibabel.ArraySequence[ndarray]
        the list of streamlines (nd.arrays) generated during tractography
    """

    # dirs shape checking
    if not (len(dirs.shape) == 4 and dirs.shape[-1] == 3):
        raise ValueError(f"dirs.shape must be of shape (M, N, P, 3), given {dirs.shape}")

    # mask shape checking
    if not (len(mask.shape) == 3 and mask.shape == dirs.shape[:3]):
        raise ValueError(f"mask.shape must be 3D and same as first 3 dimensions of dirs, given mask of shape {mask.shape} for dirs of shape {dirs.shape}")

    # valid directions (magnitude == 1 for all voxels in mask) checking
    norms = np.linalg.norm(dirs, axis=3)
    proper_norms_in_mask = int(np.sum(norms * mask))
    if proper_norms_in_mask != np.sum(mask):
        raise ValueError(f"dirs must be unit vector for all voxels within mask")
    
    # valid affine shape checking
    if not np.all(affine.shape == (4, 4)):
        raise ValueError(f"affine.shape must be (4, 4), given affine of shape {affine.shape}")

    # valid step size type checking
    if not isinstance(step_size, (float, int)):
        raise TypeError(f"step_size must be of type float or int, given {type(step_size)}")
    
    # valid max length type checking
    if not isinstance(max_length, (float, int, type(None))):
        raise TypeError(f"max_length must be of type float, int, or None, given {type(max_length)}")

    # valid max propagation angle checking
    if not isinstance(max_prop_angle, (float, int, type(None))):
        raise TypeError(f"max_prop_angle must be of type float, int, or None, given {type(max_prop_angle)}")
    
    # valid seed density checking
    if not isinstance(seed_density, int):
        raise TypeError(f"seed_density must be of type int, given {type(seed_density)}")
    
    # loading c function
    if platform.system() == "Darwin":  # macOS
        _tractography_fpath = os.path.join(os.path.dirname(__file__), "lib_tractography.dylib")
    elif platform.system() == "Windows":  # windows
        _tractography_fpath = os.path.join(os.path.dirname(__file__), "lib_tractography.dll")
    else:  # all other, .so is the meson default file type
        _tractography_fpath = os.path.join(os.path.dirname(__file__), "lib_tractography.so")
    
    # wrapping c function
    _tractography = ctypes.cdll.LoadLibrary(_tractography_fpath)
    trackStreams = _tractography.trackStreams
    trackStreams.restype = None
    trackStreams.argtypes = [ndpointer(ctypes.c_double, ndim=4), ndpointer(ctypes.c_int, ndim=3), ndpointer(ctypes.c_double, ndim=2),
                             ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
                             ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int,
                             ndpointer(ctypes.c_double, ndim=3), ndpointer(ctypes.c_int, ndim=1)]

    # creating c output arrays
    num_streamlines = 2 * (seed_density ** 3) * np.sum(mask)
    output_sls = np.zeros((num_streamlines, MAX_SL_POINTS, 3), dtype=np.double, order='C')
    output_sl_lengths = np.zeros((num_streamlines, ), dtype=np.int32, order='C')

    # TESTING
    # dirs[10, 10, 10, 2] = -99
    # mask[10, 11, 13] = 99
    # print(np.sum(mask))
    # print(MAX_SL_POINTS)
    # print(output_sls.shape)
    # tracking_mask = np.zeros([mask.shape[i] - 1 for i in range(3)], dtype=np.int32, order='C')
    # tracking_affine = np.zeros([4, 4], order='C')
    print(np.argwhere(mask == 1))
    print(affine)
    # import copy
    test_mask = np.zeros(mask.shape, dtype=np.int32, order='C')
    test_mask[18, 24, 9] = 1

    # modifying parameters for passing to c function
    if max_length is None:  # modify for c function (-1.0 is interpreted as no max length)
        c_max_length = -1.0
    if max_prop_angle is None:  # modify for c function (-1.0 is interpreted as no max prop angle)
        c_max_prop_angle = -1.0
    c_dirs = np.asarray(dirs, dtype=np.double, order='C')
    c_mask = np.asarray(mask, dtype=np.int32, order='C')  # function expects int32 data type
    c_affine = np.asarray(affine, dtype=np.double, order='C')

    # running tractography
    trackStreams(c_dirs, test_mask, c_affine,
                 dirs.shape[0], dirs.shape[1], dirs.shape[2],
                 seed_density, step_size, c_max_length, c_max_prop_angle, TractographyMethods.RUNGE_KUTTA_4,
                 output_sls, output_sl_lengths)  # testing

    # # TESTING
    # print(output_sls[0][0][0])
    # from fury import window, actor
    # scene = window.Scene()
    # scene.add(actor.contour_from_roi(mask, affine, opacity=0.5))
    # scene.add(actor.contour_from_roi(tracking_mask, tracking_affine, [0, 0, 1], opacity=0.5))
    # window.show(scene)
