# import python
MESSAGE(STATUS "Python library set: ${PYTHON_LIBRARY}")
FIND_PACKAGE(PythonInterp REQUIRED)
FIND_PACKAGE(PythonLibs REQUIRED)
FIND_PACKAGE(NumPy REQUIRED)

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

# set output directory variable for modules
IF(READDY_DEBUG_PYTHON_MODULES)
    SET(PYTHON_CPP_MODULES_OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/python/readdy/_internal")
ELSE()
    SET(PYTHON_CPP_MODULES_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/readdy/_internal")
ENDIF()

# add modules
ADD_SUBDIRECTORY(src/cxx)

# install via setup.py
FIND_PROGRAM(PYTHON "python")
if (PYTHON)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/src/python/setup.py")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    add_custom_command(OUTPUT ${OUTPUT}
            COMMAND ${PYTHON} ${SETUP_PY} build
            COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
            DEPENDS api prototyping)

    add_custom_target(target ALL DEPENDS ${OUTPUT})

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
ELSE()
    MESSAGE(FATAL_ERROR "No python executable was found!")
endif()