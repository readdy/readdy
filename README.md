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
This builds ReaDDy2 with python-bindings for python3 by default. If you want
to build for python2, set the BUILD_PY3 variable via the commandline.

	$ cmake .. -DBUILD_PY3:BOOL=OFF
You might have a certain python distribution that you want to build ReaDDy2 
for. In this case you have to specify those environment variables as well. For
example, on Ubuntu 14.04 with python2, installed via apt-get, the 
configuration could look as follows

	$ PYTHON_EXECUTABLE=/usr/bin/python2.7 PYTHON_LIBRARY=/usr/lib/python2.7/config-x86_64-linux-gnu/libpython2.7.so PYTHON_INCLUDE=/usr/include/python2.7 cmake .. -DBUILD_PY3:BOOL=OFF
The environment variables PYTHON_EXECUTABLE, PYTHON_LIBRARY and 
PYTHON_INCLUDE have to be changed according to your python distribution.

If you want to build ReaDDy2 without python-bindings to use it as a C++ 
library only, set the READDY_BUILD_PYTHON_WRAPPER variable

	$ cmake .. -DREADDY_BUILD_PYTHON_WRAPPER:BOOL=OFF