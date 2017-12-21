#####################################################################
# Copyright (c) 2016 Computational Molecular Biology Group,         #
#                    Freie Universitaet Berlin (GER)                #
#                                                                   #
# This file is part of ReaDDy.                                      #
#                                                                   #
# ReaDDy is free software: you can redistribute it and/or modify    #
# it under the terms of the GNU Lesser General Public License as    #
# published by the Free Software Foundation, either version 3 of    #
# the License, or (at your option) any later version.               #
#                                                                   #
# This program is distributed in the hope that it will be useful,   #
# but WITHOUT ANY WARRANTY; without even the implied warranty of    #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     #
# GNU Lesser General Public License for more details.               #
#                                                                   #
# You should have received a copy of the GNU Lesser General         #
# Public License along with this program. If not, see               #
# <http://www.gnu.org/licenses/>.                                   #
#####################################################################


PROJECT(runUnitTests_cpu)

FIND_PACKAGE(Threads REQUIRED)

LIST(APPEND READDY_CPU_TEST_SOURCES TestMain.cpp)
LIST(APPEND READDY_CPU_TEST_SOURCES TestKernel.cpp)
LIST(APPEND READDY_CPU_TEST_SOURCES TestNeighborList.cpp)
LIST(APPEND READDY_CPU_TEST_SOURCES TestNeighborListIterator.cpp)
LIST(APPEND READDY_CPU_TEST_SOURCES TestParallelizedGillespie.cpp)
LIST(APPEND READDY_CPU_TEST_SOURCES TestReactions.cpp)
LIST(APPEND READDY_CPU_TEST_SOURCES TestCellLinkedList.cpp)

ADD_EXECUTABLE(${PROJECT_NAME} ${READDY_CPU_TEST_SOURCES} ${TESTING_INCLUDE_DIR})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS} ${TESTING_INCLUDE_DIR} ${CPU_INCLUDE_DIR} ${GOOGLETEST_INCLUDE} ${GOOGLEMOCK_INCLUDE})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} readdy readdy_kernel_cpu gtest_main Threads::Threads)
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
        LINK_FLAGS "${EXTRA_LINK_FLAGS}"
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
ADD_TEST(NAME readdy-tests-cpu COMMAND ${PROJECT_NAME})

IF (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
    INSTALL(TARGETS ${PROJECT_NAME}
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
ENDIF (READDY_INSTALL_UNIT_TEST_EXECUTABLE)