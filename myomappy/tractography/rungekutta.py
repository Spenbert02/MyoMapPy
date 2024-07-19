"""
Runge Kutta Tractography functions
"""


import numpy as np
import os.path
from dipy.direction.peaks import PeaksAndMetrics
from dipy.tracking.utils import seeds_from_mask
from fury import window, actor


def rk4Track(pam: PeaksAndMetrics,
             mask: np.ndarray,
             affine: np.ndarray,
             step_size,
             max_length=None,
             max_prop_angle=None,
             seed_density=1):
    """generate and return an ArraySequence of streamlines (nd.arrays) created using 4th order
    Runge-Kutta tracking

    Parameters
    ----------
    peaks : dipy.direction.peaks.PeaksAndMetrics instance
        the tensor ODF PeaksAndMetrics object, used to get the primary diffusion direction for tracking
    mask : _type_
        _description_
    affine : _type_
        _description_
    step_size : _type_
        _description_
    max_length : _type_, optional
        _description_, by default None
    max_prop_angle : _type_, optional
        _description_, by default None
    seed_density : int, optional
        _description_, by default 1
    """

    # generate seeds
    seeds = seeds_from_mask(mask, affine, seed_density)
    centers = []
    init_dirs = []
    for i in range(seeds.shape[0]):
        x1 = np.array(seeds[i, :], dtype=np.double)
        x1_aff = np.array([[x1[0]], [x1[1]], [x1[2]], [1]])
        x0_aff = affine @ x1_aff
        x0 = np.array([x0_aff[0, 0], x0_aff[1, 0], x0_aff[2, 0]])
        v_x0 = np.zeros((3,), dtype=np.double)
        res = pam.get_direction(x1, v_x0)
        if res == 0:
            centers.append(x0)
            init_dirs.append(v_x0)

    return np.array(centers), np.array(init_dirs)


    # x1 = np.array([
    #     [20],
    #     [24],
    #     [9]
    # ])
    # print("in mask: ", mask[x1[0, 0], x1[1, 0], x1[2, 0]])

    # x0 = np.linalg.inv(affine) @ np.vstack([x1, [[1]]])
    # x0 = x0.reshape(4)[:3]
    # print("x0: ", x0)

    # v_x_0 = np.zeros((3,), dtype=np.double, order='C')
    # res = pam.get_direction(x0, v_x_0)
    # print("result: ", res)
    # print("v_x_0: ", v_x_0)





# --------------------------------
# The following code uses a c implementation which is not functional (Segmentation fault issues).
# --------------------------------

# import ctypes
# import platform
# from myomappy.tractography import TractographyMethods, MAX_SL_POINTS
# from myomappy.tractography._utils import unpackCArray
# from numpy.ctypeslib import ndpointer

#
# USE_C_MODULE = False  # whether to use the c module. can't get around segmentation faults, so
#                       # the python implementation will have to do for now

# def rk4Track(dirs, mask, affine, step_size, max_length=None, max_prop_angle=None, seed_density=1):
#     """from a direction field and mask, use 4th order Runge-Kutta tractography to
#     create and return a set of tracts

#     Parameters
#     ----------
#     dirs : ndarray
#         (M, N, P, 3) shape nd.array, the per-voxel direction field to track
#     mask : ndarray
#         (M, N, P) shape nd.array, the binary voxel mask to track within
#     affine : ndarray
#         (4, 4) shape nd.array, the voxel-to-world homogenous transformation matrix. must have isotropic voxel size
#     step_size : float, int
#         the step size to use during tracking
#     max_length : float, int, None
#         the (inclusive) maximum tract length. Upon reaching this length, tracts are terminated. If
#         None, there is no limit on tract length
#     max_prop_angle : float, int, None
#         the maximum tract propagation angle (degrees). If this angle is exceeded during tracking, the
#         tract is terminated. If None, propagation angle is not a termination criteria.
#     seed_density: int
#         the seed density to use for tractography. e.g., a seed density of 2 uses 2x2x2 seed points
#         per voxel.

#     Returns
#     -------
#     nibabel.ArraySequence[ndarray]
#         the list of streamlines (nd.arrays) generated during tractography
#     """

#     # dirs shape checking
#     if not (len(dirs.shape) == 4 and dirs.shape[-1] == 3):
#         raise ValueError(f"dirs.shape must be of shape (M, N, P, 3), given {dirs.shape}")

#     # mask shape checking
#     if not (len(mask.shape) == 3 and mask.shape == dirs.shape[:3]):
#         raise ValueError(f"mask.shape must be 3D and same as first 3 dimensions of dirs, given mask of shape {mask.shape} for dirs of shape {dirs.shape}")

#     # valid directions (magnitude == 1 for all voxels in mask) checking
#     norms = np.linalg.norm(dirs, axis=3)
#     proper_norms_in_mask = int(np.sum(norms * mask))
#     if proper_norms_in_mask != np.sum(mask):
#         raise ValueError(f"dirs must be unit vector for all voxels within mask")
    
#     # valid affine shape checking
#     if not np.all(affine.shape == (4, 4)):
#         raise ValueError(f"affine.shape must be (4, 4), given affine of shape {affine.shape}")

#     # TODO: here, check for isotropic voxels

#     # valid step size type checking
#     if not isinstance(step_size, (float, int)):
#         raise TypeError(f"step_size must be of type float or int, given {type(step_size)}")
    
#     # valid max length type checking
#     if not isinstance(max_length, (float, int, type(None))):
#         raise TypeError(f"max_length must be of type float, int, or None, given {type(max_length)}")

#     # valid max propagation angle checking
#     if not isinstance(max_prop_angle, (float, int, type(None))):
#         raise TypeError(f"max_prop_angle must be of type float, int, or None, given {type(max_prop_angle)}")
    
#     # valid seed density checking
#     if not isinstance(seed_density, int):
#         raise TypeError(f"seed_density must be of type int, given {type(seed_density)}")
    
#     # TODO: ensure seed density is greater or equal to 1
    
#     if USE_C_MODULE:  # this currently (7/19/2024) is not working. can't get around the fucking segmentation fault
#         # loading c function
#         if platform.system() == "Darwin":  # macOS
#             _tractography_fpath = os.path.join(os.path.dirname(__file__), "lib_tractography.dylib")
#         elif platform.system() == "Windows":  # windows
#             _tractography_fpath = os.path.join(os.path.dirname(__file__), "lib_tractography.dll")
#         else:  # all other, .so is the meson default file type
#             _tractography_fpath = os.path.join(os.path.dirname(__file__), "lib_tractography.so")
        
#         # loading c function
#         _tractography = ctypes.cdll.LoadLibrary(_tractography_fpath)
#         trackStreams = _tractography.trackStreams
#         trackStreams.restype = None
#         trackStreams.argtypes = [ndpointer(ctypes.c_double, ndim=4), ndpointer(ctypes.c_int, ndim=3), ndpointer(ctypes.c_double, ndim=2),
#                                 ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
#                                 ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int,
#                                 ndpointer(ctypes.c_double, ndim=3), ndpointer(ctypes.c_int, ndim=1)]

#         # creating c output arrays
#         num_streamlines = 2 * (seed_density ** 3) * np.sum(mask)
#         output_sls = np.zeros((num_streamlines, MAX_SL_POINTS, 3), dtype=np.double, order='C')
#         output_sl_lengths = np.zeros((num_streamlines, ), dtype=np.int32, order='C')

#         # modifying parameters for passing to c function
#         if max_length is None:  # modify for c function (-1.0 is interpreted as no max length)
#             c_max_length = -1.0
#         if max_prop_angle is None:  # modify for c function (-1.0 is interpreted as no max prop angle)
#             c_max_prop_angle = -1.0
#         c_dirs = np.asarray(dirs, dtype=np.double, order='C')
#         c_mask = np.asarray(mask, dtype=np.int32, order='C')  # function expects int32 data type
#         c_affine = np.asarray(affine, dtype=np.double, order='C')

#         # running tractography
#         trackStreams(c_dirs, c_mask, c_affine,
#                     dirs.shape[0], dirs.shape[1], dirs.shape[2],
#                     seed_density, step_size, c_max_length, c_max_prop_angle, TractographyMethods.EULER,  # NOTE: for testing, using euler tracking
#                     output_sls, output_sl_lengths)
        
#         ret_streams = unpackCArray(output_sls, output_sl_lengths)
#         return(ret_streams)
