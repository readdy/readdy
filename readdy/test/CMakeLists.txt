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

project(runUnitTests)

add_executable(${PROJECT_NAME}
        TestMain.cpp TestAlgorithms.cpp TestCompartments.cpp TestContext.cpp
        TestDetailedBalance.cpp TestIndex.cpp TestIndexPersistentVector.cpp TestIntegration.cpp
        TestMatrix33.cpp TestObservables.cpp TestPlugins.cpp TestPotentials.cpp TestReactions.cpp
        TestSignals.cpp TestSimulationLoop.cpp TestStateModel.cpp TestTopologies.cpp TestTopologyGraphs.cpp
        TestTopologyReactions.cpp TestTopologyReactionsExternal.cpp TestVec3.cpp
        ${TESTING_INCLUDE_DIR})

target_include_directories(
        ${PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS} ${TESTING_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} PUBLIC readdy Catch2::Catch2)

set_target_properties(${PROJECT_NAME} PROPERTIES
        LINK_FLAGS "${EXTRA_LINK_FLAGS}"
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
catch_discover_tests(${PROJECT_NAME})

if (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
    install(TARGETS ${PROJECT_NAME}
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
endif (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
