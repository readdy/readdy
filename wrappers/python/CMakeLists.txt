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
            DEPENDS api prototyping common)

    add_custom_target(target ALL DEPENDS ${OUTPUT})

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
ELSE()
    MESSAGE(FATAL_ERROR "No python executable was found!")
endif()