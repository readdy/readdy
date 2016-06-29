IF (NOT NUMPY_INCLUDE_DIR)
    FIND_PACKAGE(PythonInterp REQUIRED)
    IF (PYTHON_EXECUTABLE)
        EXECUTE_PROCESS(
                COMMAND ${PYTHON_EXECUTABLE} "${READDY_GLOBAL_DIR}/libraries/boost/numpy_include_dir.py"
                OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
                RESULT_VARIABLE Result
        )
        MESSAGE(STATUS "Found numpy include dir: ${NUMPY_INCLUDE_DIR}, result=${Result}")
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