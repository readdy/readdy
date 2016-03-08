# Find pthread library
FIND_PACKAGE(Threads REQUIRED)

LINK_DIRECTORIES(${GOOGLETEST_LIBRARY})
LINK_LIBRARIES(gtest gtest_main ${CMAKE_THREAD_LIBS_INIT} ReaDDy2)


FILE(GLOB READDY_TEST_SOURCES *.cpp)
ADD_EXECUTABLE(runUnitTests ${READDY_TEST_SOURCES} ${GOOGLETEST_INCLUDE})

INCLUDE_DIRECTORIES(BEFORE ${ReaDDy2_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(runUnitTests ReaDDy2 ${gtest})

ADD_TEST(NAME readdy-tests COMMAND runUnitTests)


#
# memory tests
#
IF(READDY_RUN_TESTS_MEMORY_CHECK)
    FIND_PACKAGE(Valgrind QUIET)
    IF (NOT VALGRIND_FOUND)
        MESSAGE(STATUS "Valgrind not found, therefore no memory checks.")
    ELSE (NOT VALGRIND_FOUND)
        MESSAGE(STATUS "Valgrind found, memory checks enabled.")
        SET(VALGRIND_OPTIONS --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=30 --leak-check=summary --leak-resolution=high)

        FUNCTION(ADD_MEMCHECK_TEST name binary)
            SET(memcheck_command "valgrind ${VALGRIND_OPTIONS}")
            SEPARATE_ARGUMENTS(memcheck_command)
            ADD_TEST(memcheck_${name} ${memcheck_command} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${binary}" ${ARGN})
        ENDFUNCTION(ADD_MEMCHECK_TEST)

        # can be called by "ctest -V -R memcheck_readdy-tests"
        ADD_MEMCHECK_TEST(readdy-tests runUnitTests)
    ENDIF (NOT VALGRIND_FOUND)
ENDIF()