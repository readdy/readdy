####################################################################
# Copyright © 2018 Computational Molecular Biology Group,          #
#                  Freie Universität Berlin (GER)                  #
#                                                                  #
# Redistribution and use in source and binary forms, with or       #
# without modification, are permitted provided that the            #
# following conditions are met:                                    #
#  1. Redistributions of source code must retain the above         #
#     copyright notice, this list of conditions and the            #
#     following disclaimer.                                        #
#  2. Redistributions in binary form must reproduce the above      #
#     copyright notice, this list of conditions and the following  #
#     disclaimer in the documentation and/or other materials       #
#     provided with the distribution.                              #
#  3. Neither the name of the copyright holder nor the names of    #
#     its contributors may be used to endorse or promote products  #
#     derived from this software without specific                  #
#     prior written permission.                                    #
#                                                                  #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           #
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      #
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         #
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         #
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            #
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     #
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         #
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER #
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      #
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       #
####################################################################

project(runUnitTests)

add_executable(${PROJECT_NAME}
        TestMain.cpp TestAlgorithms.cpp TestCompartments.cpp TestContext.cpp
        TestDetailedBalance.cpp TestIndex.cpp TestIndexPersistentVector.cpp TestIntegration.cpp
        TestMatrix33.cpp TestObservables.cpp TestPlugins.cpp TestPotentials.cpp TestReactions.cpp
        TestSignals.cpp TestSimulationLoop.cpp TestStateModel.cpp TestTopologies.cpp TestTopologyGraphs.cpp
        TestTopologyReactions.cpp TestTopologyReactionsExternal.cpp TestVec3.cpp TestBreakingBonds.cpp IntegrationTests.cpp
        TestGeometries.cpp Graph.cpp IndexPersistentVector.cpp Vertex.cpp H5RD.cpp
        ${TESTING_INCLUDE_DIR})

if (READDY_BUILD_MPI_KERNEL)
    find_package(MPI REQUIRED)
    # MPI specific flags for target library
    target_include_directories(${PROJECT_NAME} PUBLIC ${MPI_CXX_INCLUDE_DIRS})
    target_compile_options(${PROJECT_NAME} PUBLIC ${MPI_CXX_COMPILE_OPTIONS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
endif()

target_include_directories(
        ${PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS} ${TESTING_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} PUBLIC readdy Catch2::Catch2)

set_target_properties(${PROJECT_NAME} PROPERTIES
        LINK_FLAGS "${EXTRA_LINK_FLAGS}"
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")

if (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
    install(TARGETS ${PROJECT_NAME}
            RUNTIME DESTINATION ${INSTALL_DIR}
            LIBRARY DESTINATION ${INSTALL_DIR}
            ARCHIVE DESTINATION ${INSTALL_DIR})
endif (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
