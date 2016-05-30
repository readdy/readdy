### Status
| Travis | Appveyor |
| --- | --- |
|[![Build Status](https://travis-ci.org/readdy/readdy.svg?branch=master)](https://travis-ci.org/readdy/readdy) | [![Build status](https://ci.appveyor.com/api/projects/status/ve6rhy2fs2fnyjqi?svg=true)](https://ci.appveyor.com/project/clonker/readdy2) |

### Dependencies
- HDF5
- boost
- cmake
- *optional*: python (2 or 3)

### Project structure
```
ReaDDy2/
|   README.md
|   ...
|
|___bin/
|   |   (binary files)
|   |
|   |___debug/
|   |   |   (debug binaries)
|   |
|   |___release/
|   |   |   (release binaries)
|
|___kernels/
|   |
|___|___include/
|   |   |   *.h (low level api || does this work? investigate on plugin structures in c++)
|   |
|   |___singlecpu/
|   |   |   ...
|   |
|   |___cuda/
|   |   |   ...
|
|___include/
|   |   *.h (public c++ api)
|
|___readdy2/
|   |   (readdy2 source code + private headers)
|   |
|___|___io/
|___|___|___tests/
|   |   |   |   (io specific testing)
|   |
|___|___adapter/
|___|___|___tests/
|   |   |   |   (adapter specific testing)
|   |   
|___|___tests/
|   |   |   (integration tests)
|
|___wrappers/
|___|___python/
|   |   |   *.py
|   |   |
|___|___|___tests/
|   |   |   |   *.py
|
|___libraries/
|   |   (directory structure with libraries)

```
### Building
Build in-source by typing:

	$ mkdir build
	$ cd build
	$ cmake .. 
	$ make
This builds ReaDDy2 with python-bindings. If you have multiple python distributions
you might run into problems. In this case you have to specify the
environment variables of the desired python distribution as well. For
example, on Ubuntu 14.04 with python2, installed via apt-get, the 
configuration could look as follows

	$ PYTHON_EXECUTABLE=/usr/bin/python2.7 \
		PYTHON_LIBRARY=/usr/lib/python2.7/config-x86_64-linux-gnu/libpython2.7.so \
		PYTHON_INCLUDE=/usr/include/python2.7 \
		cmake .. -DREADDY_BUILD_PYTHON_WRAPPER:BOOL=ON
The environment variables PYTHON_EXECUTABLE, PYTHON_LIBRARY and 
PYTHON_INCLUDE have to be changed according to your python distribution.

If you want to build ReaDDy2 without python-bindings to use it as a C++ 
library only, set the READDY_BUILD_PYTHON_WRAPPER OFF variable

	$ cmake .. -DREADDY_BUILD_PYTHON_WRAPPER:BOOL=OFF
	
Deciding if shared or static: Set CMake parameter BUILD_SHARED_LIBS.