####################################################################
# Copyright © 2020 Computational Molecular Biology Group,          #
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

project(run_readdy_scenarios)

#include("${READDY_GLOBAL_DIR}/cmake/sources/kernels/cpu.cmake")

list(APPEND SCENARIO_SOURCES "main.cpp")
list(APPEND SCENARIO_LIBS "readdy")
#list(APPEND SCENARIO_LIBS "readdy_kernel_cpu")
list(APPEND SCENARIO_INCLUDE_DIRS "${READDY_INCLUDE_DIRS}")
list(APPEND SCENARIO_INCLUDE_DIRS "${CPU_INCLUDE_DIR}")

#if(READDY_BUILD_MPI_KERNEL)
#    # ok, will link to mpi indirectly through readdy_kernel_mpi
#else()
#    # the MPISession object makes some MPI_* calls that need to be mocked if not running MPI
#    list(APPEND SCENARIO_SOURCES "${CMAKE_SOURCE_DIR}/readdy_testing/mpimock/mpi.cpp")
#    list(APPEND SCENARIO_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/readdy_testing/mpimock")
#endif()

add_executable(${PROJECT_NAME} ${SCENARIO_SOURCES})
message(STATUS ${SCENARIO_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${SCENARIO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${SCENARIO_LIBS})
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")

if (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
    install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
endif (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
