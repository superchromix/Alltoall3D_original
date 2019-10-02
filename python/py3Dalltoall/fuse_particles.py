"""
    Python binding for 3Dalltoall.
    See https://github.com/berndrieger/3Dalltoall

    The binding is based on ctypes.
    See https://docs.python.org/3.5/library/ctypes.html, http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
"""

import os
from ctypes import cdll, POINTER, c_int32, c_double
import numpy as np

# define library loader (actual loading is lazy)
package_dir = os.path.dirname(os.path.realpath(__file__))

if os.name == 'nt':
    lib_path = os.path.join(package_dir, 'fuse_particles_3d.dll') # library name on Windows
elif os.name == 'posix':
    lib_path = os.path.join(package_dir, 'libfuse_particles_3d.so') # library name on Unix
else:
    raise RuntimeError('OS {} not supported by py3Dalltoall.'.format(os.name))

lib = cdll.LoadLibrary(lib_path)

# fuse_particles_3d function in the dll
func = lib.fuse_particles_3d
func.restype = c_int32
func.argtypes = [
    POINTER(c_double),  # transformed_coordinates_x
    POINTER(c_double),  # transformed_coordinates_y
    POINTER(c_double),  # transformed_coordinates_z
    POINTER(c_double),  # transformation_parameters
    c_int32,  # number_particles
    POINTER(c_int32),  # localizations_per_particle
    POINTER(c_double),  # coordinates_x
    POINTER(c_double),  # coordinates_y
    POINTER(c_double),  # coordinates_z
    POINTER(c_double),  # precision_xy
    POINTER(c_double),  # precision_z
    c_double,  # gauss_transform_scale
    POINTER(c_int32),  # channel_ids
    c_int32,  # averaging_channel_id
    c_int32,  # number_iterations_all2all
    c_int32,  # number_iterations_one2all
    c_int32,  # symmetry_order
    c_double  # outlier_threshold
]


def fuse_particles_3d(localizations_per_particle, coordinates_x, coordinates_y, coordinates_z, precision_xy, precision_z,
                      gauss_transform_scale, channel_ids=None, averaging_channel_id=0, number_iterations_all2all=1, number_iterations_one2all=10,
                      symmetry_order=0, outlier_threshold=1):
    """

    :return:
    """

    # checks
    number_particles = localizations_per_particle.size
    number_localizations = np.sum(localizations_per_particle)
    d = (number_localizations,)

    if channel_ids is None:
        channel_ids = np.zeros(d, dtype=np.int32)

    if any([not x.flags.c_contiguous for x in [coordinates_x, coordinates_y, coordinates_z, precision_xy, precision_z, channel_ids]]):
        raise RuntimeError('Memory layout of data arrays mismatch')

    # pre-allocate output variables
    d = (number_localizations * (number_iterations_one2all + 1),)
    transformed_coordinates_x = np.zeros(d, dtype=np.double)
    transformed_coordinates_y = np.zeros(d, dtype=np.double)
    transformed_coordinates_z = np.zeros(d, dtype=np.double)
    transformation_parameters = np.zeros((16 * number_particles * (number_iterations_one2all + 1), ), dtype=np.double)

    # call into the library
    status = func(
        transformed_coordinates_x.ctypes.data_as(func.argtypes[0]),
        transformed_coordinates_y.ctypes.data_as(func.argtypes[1]),
        transformed_coordinates_z.ctypes.data_as(func.argtypes[2]),
        transformation_parameters.ctypes.data_as(func.argtypes[3]),
        func.argtypes[4](number_particles),
        localizations_per_particle.ctypes.data_as(func.argtypes[5]),
        coordinates_x.ctypes.data_as(func.argtypes[6]),
        coordinates_y.ctypes.data_as(func.argtypes[7]),
        coordinates_z.ctypes.data_as(func.argtypes[8]),
        precision_xy.ctypes.data_as(func.argtypes[9]),
        precision_z.ctypes.data_as(func.argtypes[10]),
        func.argtypes[11](gauss_transform_scale),
        channel_ids.ctypes.data_as(func.argtypes[12]),
        func.argtypes[13](averaging_channel_id),
        func.argtypes[14](number_iterations_all2all),
        func.argtypes[15](number_iterations_one2all),
        func.argtypes[16](symmetry_order),
        func.argtypes[17](outlier_threshold)
    )

    # check status
    if status != 0:
        raise RuntimeError('status = {}'.format(status))

    # return output values
    return transformed_coordinates_x, transformed_coordinates_y, transformed_coordinates_z, transformation_parameters
