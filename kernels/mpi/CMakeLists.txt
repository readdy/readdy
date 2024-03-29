####################################################################
# Copyright © 2019 Computational Molecular Biology Group,          #
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

project(readdy_kernel_mpi C CXX)

find_package(MPI REQUIRED)

# sources and headers
include("${READDY_GLOBAL_DIR}/cmake/sources/kernels/mpi.cmake")

# create library
add_library(${PROJECT_NAME} SHARED ${MPI_SOURCES} ${READDY_INCLUDE_DIRS} ${MPI_INCLUDE_DIR})

# MPI specific flags for target library
target_include_directories(${PROJECT_NAME} PUBLIC ${MPI_CXX_INCLUDE_DIRS})
target_compile_options(${PROJECT_NAME} PUBLIC ${MPI_CXX_COMPILE_OPTIONS})
target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})

# includes
target_include_directories(${PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS} ${MPI_INCLUDE_DIR})

# link
target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_DL_LIBS} readdy)

set_target_properties(${PROJECT_NAME} PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        LINK_FLAGS "${EXTRA_LINK_FLAGS}"
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
        ARCHIVE_OUTPUT_DIRECTORY "${READDY_PLUGIN_OUTPUT_DIRECTORY}"
        LIBRARY_OUTPUT_DIRECTORY "${READDY_PLUGIN_OUTPUT_DIRECTORY}"
        LINKER_LANGUAGE CXX)

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION ${INSTALL_DIR}/readdy_plugins
        LIBRARY DESTINATION ${INSTALL_DIR}/readdy_plugins
        ARCHIVE DESTINATION ${INSTALL_DIR}/readdy_plugins
        )

if (READDY_CREATE_TEST_TARGET)
    add_subdirectory(test)
endif()