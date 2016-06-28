#  NUMPY_INCLUDE_DIR, where to find numpy/arrayobject.h, etc.
#  NUMPY_FOUND, If false, do not try to use numpy headers.
IF (NOT NUMPY_INCLUDE_DIR)
    FIND_PACKAGE(PythonInterp REQUIRED)
    IF (PYTHON_EXECUTABLE)
        EXECUTE_PROCESS(
                COMMAND ${PYTHON_EXECUTABLE} "-c 'import numpy; print numpy.get_include()'"
                OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
                RESULT_VARIABLE NUMPY_NOT_FOUND
        )
        IF(NOT Result EQUAL "0")
            MESSAGE(WARNING "No NumPy installed for interpreter ${PYTHON_EXECUTABLE}.")
            SET(NUMPY_FOUND FALSE)
        ELSE()
            SET(NUMPY_FOUND TRUE)
            SET (NUMPY_INCLUDE_DIR ${NUMPY_INCLUDE_DIR} CACHE PATH "Numpy include path")
        ENDIF()

    ELSE()
        MESSAGE(STATUS "FindNumPy: Python executable not found")
    ENDIF()
endif (NOT NUMPY_INCLUDE_DIR)