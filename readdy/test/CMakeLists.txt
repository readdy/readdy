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

ADD_EXECUTABLE(runUnitTests
        TestMain.cpp TestAggregators.cpp TestAlgorithms.cpp TestCompartments.cpp TestContext.cpp
        TestDetailedBalance.cpp TestIndex.cpp TestIndexPersistentVector.cpp TestIntegration.cpp
        TestMatrix33.cpp TestObservables.cpp TestPlugins.cpp TestPotentials.cpp TestReactions.cpp
        TestSignals.cpp TestSimulationSchemes.cpp TestStateModel.cpp TestTopologies.cpp TestTopologyGraphs.cpp
        TestTopologyReactions.cpp TestTopologyReactionsExternal.cpp TestVec3.cpp
        ${TESTING_INCLUDE_DIR})

TARGET_INCLUDE_DIRECTORIES(
        runUnitTests PUBLIC ${READDY_INCLUDE_DIRS} ${TESTING_INCLUDE_DIR}
)

TARGET_COMPILE_DEFINITIONS(runUnitTests PRIVATE READDY_KERNELS_TO_TEST="${READDY_KERNELS_TO_TEST}")

TARGET_LINK_LIBRARIES(runUnitTests PUBLIC readdy Catch2::Catch2)

SET_TARGET_PROPERTIES(runUnitTests PROPERTIES
        LINK_FLAGS "${EXTRA_LINK_FLAGS}"
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
catch_discover_tests(runUnitTests)

IF (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
    INSTALL(TARGETS runUnitTests
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
ENDIF (READDY_INSTALL_UNIT_TEST_EXECUTABLE)

# CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/simple_trajectory.h5 ${CMAKE_CURRENT_BINARY_DIR}/simple_trajectory.h5)
