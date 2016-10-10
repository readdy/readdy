setlocal enableextensions

%PYTHON% -c "from __future__ import print_function; import distutils.sysconfig; print(distutils.sysconfig.get_python_inc(True))" > python_include_dir.txt
set /p PYTHON_INCLUDE_DIR=<python_include_dir.txt

set CMAKE_FLAGS=-DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX%
set CMAKE_FLAGS=-DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX%
set CMAKE_FLAGS=%CMAKE_FLAGS% -DREADDY_BUILD_SHARED_COMBINED:BOOL=ON
set CMAKE_FLAGS=%CMAKE_FLAGS% -DCMAKE_BUILD_TYPE=Release
set CMAKE_FLAGS=%CMAKE_FLAGS% -DHDF5_INCLUDE_DIR=%LIBRARY_PREFIX%\include

set HDF5_ROOT=%LIBRARY_PREFIX%

mkdir build
cd build

cmake .. -G "Visual Studio 14 2015 Win64" %CMAKE_FLAGS%
if %errorlevel% neq 0 exit /b %errorlevel%

cmake --build . --config Release --target INSTALL
