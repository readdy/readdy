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


# Find pthread library
FIND_PACKAGE(Threads REQUIRED)

LIST(APPEND READDY_TEST_SOURCES TestMain.cpp)
LIST(APPEND READDY_TEST_SOURCES TestPlugins.cpp)
LIST(APPEND READDY_TEST_SOURCES TestReactions.cpp)
LIST(APPEND READDY_TEST_SOURCES TestSimulation.cpp)
LIST(APPEND READDY_TEST_SOURCES TestStateModel.cpp)
LIST(APPEND READDY_TEST_SOURCES TestPotentials.cpp)
LIST(APPEND READDY_TEST_SOURCES TestObservables.cpp)
LIST(APPEND READDY_TEST_SOURCES TestSignals.cpp)
LIST(APPEND READDY_TEST_SOURCES TestKernelContext.cpp)
LIST(APPEND READDY_TEST_SOURCES TestSimulationSchemes.cpp)
LIST(APPEND READDY_TEST_SOURCES TestCompartments.cpp)
LIST(APPEND READDY_TEST_SOURCES TestVec3.cpp)
LIST(APPEND READDY_TEST_SOURCES TestAggregators.cpp)

IF(INCLUDE_PERFORMANCE_TESTS)
    LIST(APPEND READDY_TEST_SOURCES TestPerformance.cpp)
ENDIF(INCLUDE_PERFORMANCE_TESTS)

ADD_EXECUTABLE(runUnitTests ${READDY_TEST_SOURCES} ${TESTING_INCLUDE_DIR})

TARGET_INCLUDE_DIRECTORIES(
        runUnitTests PUBLIC ${READDY_INCLUDE_DIRS} ${TESTING_INCLUDE_DIR} "${gtest_SOURCE_DIR}/include" "${gmock_SOURCE_DIR}/include"
)

TARGET_COMPILE_DEFINITIONS(runUnitTests PRIVATE READDY_KERNELS_TO_TEST="${READDY_KERNELS_TO_TEST}")

TARGET_LINK_LIBRARIES(runUnitTests PUBLIC gmock_main readdy Threads::Threads)

SET_TARGET_PROPERTIES(runUnitTests PROPERTIES
        LINK_FLAGS "${EXTRA_LINK_FLAGS}"
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")

IF (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
    INSTALL(TARGETS runUnitTests
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
ENDIF (READDY_INSTALL_UNIT_TEST_EXECUTABLE)