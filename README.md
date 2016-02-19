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
	$ cmake .. -DBUILD_PY3:BOOL=ON
	$ make
To specify if building for python3, hand in the flag BUILD_PY3 as commandline argument.