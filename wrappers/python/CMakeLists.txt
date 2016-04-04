# add READDY_WITH_PYTHON define
ADD_DEFINITIONS(-DREADDY_WITH_PYTHON)

# import python
MESSAGE(STATUS "Python library set: ${PYTHON_LIBRARY}")
FIND_PACKAGE(PythonInterp REQUIRED)
FIND_PACKAGE(PythonLibs REQUIRED)

IF (${PYTHON_VERSION_STRING} GREATER 3.0)
    MESSAGE(STATUS "--------------------- Building with python 3 ---------------------")
ELSE ()
    MESSAGE(STATUS "--------------------- Building with python 2 ---------------------")
ENDIF ()

FIND_PACKAGE(Boost COMPONENTS system thread REQUIRED)
# macro that takes a list of potential library names
MACRO(boost_python_libraries libnames)
    FOREACH (_boost_py_lib ${libnames})
        MESSAGE(STATUS "Looking for boost lib: ${_boost_py_lib}")
        FIND_PACKAGE(Boost QUIET COMPONENTS ${_boost_py_lib})
        STRING(TOUPPER ${_boost_py_lib} boost_py_lib_name)
        IF (Boost_${boost_py_lib_name}_FOUND)
            SET(Boost_PYTHON_LIBRARIES ${Boost_${boost_py_lib_name}_LIBRARIES})
            BREAK()
        ENDIF ()
    ENDFOREACH ()
ENDMACRO(boost_python_libraries)

# look for boost python
IF (${PYTHON_VERSION_STRING} GREATER 3.0)
    # find boost-python3 lib
    boost_python_libraries("python3;python-py35;python-py34;python")
ELSE ()
    # find boost-python2 lib
    boost_python_libraries("python2;python-py27;python")
ENDIF ()
# if no boost-python, exit
IF (NOT Boost_PYTHON_LIBRARIES)
    MESSAGE(FATAL_ERROR "Did not find boost python libraries, exit!")
ENDIF ()

MESSAGE(STATUS "Python libraries: ${PYTHON_LIBRARIES}")
MESSAGE(STATUS "Boost libraries: ${Boost_LIBRARIES}")
MESSAGE(STATUS "Python version: ${PYTHON_VERSION_STRING}")
MESSAGE(STATUS "Boost python library: ${Boost_PYTHON_LIBRARIES}")

PYTHON_ADD_MODULE(simulation PySimulation.cpp)
TARGET_INCLUDE_DIRECTORIES(simulation PUBLIC $ENV{PYTHON_INCLUDE_DIR})
TARGET_INCLUDE_DIRECTORIES(simulation PUBLIC "${ReaDDy2_INCLUDE_DIRS};${Boost_INCLUDE_DIRS}")
TARGET_LINK_LIBRARIES(simulation ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} ReaDDy2)

SET_TARGET_PROPERTIES(simulation
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/modules"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/modules"
)