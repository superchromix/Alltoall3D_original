# 3Dalltoall

This code is created for a Windows or Linux environment. The main code is written in MATLAB and some of the compute-intensive kernels have been written in CUDA and C++. The Matlab code gets compiled into a library using
the [Matlab compiler](https://www.mathworks.com/products/compiler.html), which can be passed to anyone without
the need to acquire a Matlab license (only requires the [Matlab runtime](https://www.mathworks.com/products/compiler/matlab-runtime.html)). The compiled Matlab code is accessed from a C interface, which is used
by the included Python package.

## Requirements

The build framework requires CMake.
If the GPU shall be used, an installed [CUDA toolkit](https://developer.nvidia.com/cuda-downloads) and the [CUB library](https://nvlabs.github.io/cub/) is required.
For the Matlab compilation, Matlab must be installed (.
For the Python packaging, a Python environment with some basic Python packages (setuptools, python-wheel) must be available.

## Build instructions

### Get the sources

The Git repository uses submodules. Include them in a _git clone_ action using the _--recursive_ option.

> git clone https://github.com/berndrieger/3Dalltoall.git --recursive

In the following

- BUILD_DIRECTORY is the directory where the project will be built
- SOURCE_DIRECTORY is the root directory of the sources
- CUB_DIRECTORY is the root directory of the downloaded [CUB library](https://nvlabs.github.io/cub/) sources

### Windows

Either use the CMake GUI or use CMake from the command line. On the command line, use the following command.

> cd BUILD_DIRECTORY<br>
> cmake -G "GENERATOR_NAME" -DCUB_ROOT_DIR=CUB_DIRECTORY SOURCE_DIRECTORY

Here, GENERATOR_NAME is the name of the target compiler/IDE and platform, e.g. Visual Studio 14 2015 Win 64 for Visual Studio 2015 x64 ([CMake Generators](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html)).

Test the example by running the example_fuse_particles_3d target. 

Install the Python binding in an active Python environment by installing the created wheel.

> pip install BUILD_DIRECTORY/py3Dalltoall/dist/py3Dalltoall-1.0.0-py2.py3-none-any.whl<br>

Test the Python binding by executing an example script.

> cd SOURCE_DIRECTORY/python/examples<br>
> python three_particles.py

### Linux

Either use the CMake GUI or use CMake from the command line. On the command line, use the following command.

> cd BUILD_DIRECTORY<br>
> cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_C_COMPILER=GCC_COMPILER_VERSION -DCUB_ROOT_DIR=CUB_DIRECTORY SOURCE_DIRECTORY

where GCC_COMPILER_VERSION is the desired C compiler version, e.g. "gcc-5" (use this CMake option if the CUDA compiler requires a gcc versions different from the default gcc version on your system, e.g. CUDA 9.1 requires a gcc version below 6).

Then build the project

> make

To use the mex files in Matlab and to use the C interface library the library paths have to be adapted

> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:MATLAB_PATH/runtime/glnxa64:MATLAB_PATH/bin/glnxa64:MATLAB_PATH/sys/os/glnxa64:MATLAB_PATH/sys/opengl/lib/glnxa64:BUILD_DIRECTORY/c_interface<br>
> export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6

where MATLAB_PATH is the root path of the installed MATLAB, e.g. "/usr/local/MATLAB/R2017b".

Test the build by running the example.

> cd BUILD_DIRECTORY/c_interface/examples<br>
> ./example_fuse_particles_3d

> virtualenv test<br>
> source test/bin/activate<br>
> pip install numpy<br>
> pip install BUILD_DIRECTORY/py3Dalltoall/dist/py3Dalltoall-1.0.0-py2.py3-none-any.whl<br>


Install the Python binding in an active Python environment by installing the created wheel.

> pip install BUILD_DIRECTORY/py3Dalltoall/dist/py3Dalltoall-1.0.0-py2.py3-none-any.whl<br>

Test the Python binding

> cd SOURCE_DIRECTORY/python/examples<br>
> python three_particles.py


### Matlab interface

The Matlab interface of the project is the _fuse_particles_3d.m_ Matlab script located in the Matlab folder.

> function [transformed_coordinates_x, transformed_coordinates_y, transformed_coordinates_z, transformation_parameters]
    = fuse_particles_3d(n_particles, n_localizations_per_particle, coordinates_x, coordinates_y,
        coordinates_z, weights_xy, weights_z, channel_ids, averaging_channel_id, n_iterations_all2all,
        n_iterations_one2all, symmetry_order, outlier_threshold)


### C interface

Using the C interface does not require Matlab to be installed, however a [Matlab Runtime](https://www.mathworks.com/products/compiler/matlab-runtime.html)
must be installed.

> int fuse_particles_3d(
        double * transformed_coordinates_x,
        double * transformed_coordinates_y,
        double * transformed_coordinates_z,
        double * transformation_parameters,
        int n_particles,
        int * n_localizations_per_particle,
        double * coordinates_x,
        double * coordinates_y,
        double * coordinates_z,
        double * weights_xy,
        double * weights_z,
        int * channel_ids,
        int averaging_channel_id,
        int n_iterations_alltoall,
        int n_iterations_onetoall,
        int symmetry_order,
        double outlier_threshold);

#### Python binding to the C interface 

Using the Python binding additionally requires [NumPy](https://www.numpy.org/). During the building a Python
wheel with the _py3Dalltoall_ module is created and can be installed using pip. The building may requires an additional install of the Python packages setuptools and wheel.

> def fuse_particles_3d(localizations_per_particle, coordinates_x, coordinates_y, coordinates_z, weights_xy, weights_z,
                      channel_ids=None, averaging_channel_id=0, number_iterations_all2all=1, number_iterations_one2all=10,
                      symmetry_order=0, outlier_threshold=1):

[//]: # (See also https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)
