## ReaDDy

ReaDDy is a Python / C++ based particle reaction-diffusion simulator. Its documentation can be found [here](https://readdy.github.io/readdy_documentation/).
*Currently, it is under construction and should not be expected to be fully functional.* For a working and more feature complete software, please see the [Java implementation of ReaDDy](https://github.com/readdy/readdy_java).

### Status
| Travis | <span style="color:gray;">Appveyor</span> |
| --- | --- |
|[![Build Status](https://travis-ci.org/readdy/readdy.svg?branch=master)](https://travis-ci.org/readdy/readdy) | not yet supported |

### Dependencies
- HDF5
- cmake
- *optional*: python (2 or 3), numpy (for python bindings)
- *testing*: gtest (included by git submodule)

### Project structure
```
readdy/
|   README.md
|   ...
|
|___kernels/
|   |
|   |___cpu/
|___|___|___include/
|   |   |   (kernel includes)
|___|___|___src/
|   |   |   (kernel sources)
|___|___|___test/
|   |   |   (kernel tests)
|   |
|   |___cuda/
|   |   |   (yet to be added)
|
|___include/
|   |   *.h (core library includes)
|
|___readdy/
|   |
|___|___main/
|   |   |   (core library sources)
|___|___test/
|   |   |   (core library tests)
|
|___libraries/
|   |   (googletest, h5xx, pybind11, spdlog)
|
|___wrappers/
|___|___python/
|   |   |   (code for python api)
|___|___|___src/
|___|___|___|___cxx/
|   |   |   |   |   (c++ code for the python module)
|   |   |   |
|___|___|___|___python/
|   |   |   |   |   (python top-level api)

```
### Building

Build by conda-build:
> Install conda-build and then execute the conda recipe:
```bash
conda install conda-build
conda-build PATH_TO_READDY/tools/conda-recipe
```

Build by using CMake: This type of build is suggested if one is interested in development of the software. There are a number of CMake options that influence the type of build:

| CMake option | Default value | Description |
| --- | --- | --- |
| READDY_CREATE_TEST_TARGET:BOOL | ON | Determining if the test targets should be generated. |
| READDY_CREATE_MEMORY_CHECK_TEST_TARGET:BOOL | OFF | Determining if the test targets should be additionally called through valgrind. Requires the previous option to be ON and valgrind to be installed. |
| READDY_INSTALL_UNIT_TEST_EXECUTABLE:BOOL | OFF | Determining if the unit test executables should be installed. This is option is mainly important for the conda recipe. |
| READDY_BUILD_SHARED_COMBINED:BOOL | OFF | Determining if the core library should be built monolithically or as separated shared libraries. |
| READDY_BUILD_PYTHON_WRAPPER:BOOL | ON | Determining if the python wrapper should be built. |
| READDY_DEBUG_PYTHON_MODULES:BOOL | OFF | If this flag is set to ON, the generated python module will be placed in-source rather than in the output directory to enable faster development. |
| READDY_DEBUG_CONDA_ROOT_DIR:PATH | "" | This option is to be used in conjunction with the previous option and only has effect if it is set to ON. It should point to the conda environment which is used for development and then effects the output directory of the binary files such that they get compiled directly into the respective environment. |
| READDY_GENERATE_DOCUMENTATION_TARGET:BOOL | OFF | Determines if the documentation target should be generated or not, which, if generated, can be called by "make doc". |
| READDY_GENERATE_DOCUMENTATION_TARGET_ONLY:BOOL | OFF | This option has the same effect as the previous option, just that it does not need any dependencies other than doxygen to be fulfilled and generates the documentation target exclusively. |
| READDY_LOG_CMAKE_CONFIGURATION:BOOL | OFF | This option determines if the status of relevant cmake cache variables should be logged at configuration time or not. |
| READDY_KERNELS_TO_TEST:STRING | "SingleCPU,CPU" | Comma separated list of kernels against which the core library should be tested within the test targets. |
| *advanced*: INCLUDE_PERFORMANCE_TESTS:BOOL | OFF | Flag indicating if the performance tests should be part of the unit test target or not. |

After having configured the cmake cache variables, one can invoke cmake and make and compile the project.

Altogether, a shell script invoking cmake with modified parameters in an environment with multiple python versions could look like [this](tools/dev/configure.sh).
