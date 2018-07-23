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

PROJECT(readdy_kernel_cpu C CXX)

# sources and headers
INCLUDE("${READDY_GLOBAL_DIR}/cmake/sources/kernels/cpu.cmake")

# create library
ADD_LIBRARY(${PROJECT_NAME} SHARED ${CPU_SOURCES} ${READDY_INCLUDE_DIRS} ${CPU_INCLUDE_DIR})

# includes
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS} ${CPU_INCLUDE_DIR})

# link
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CMAKE_DL_LIBS} readdy)

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        LINK_FLAGS "${EXTRA_LINK_FLAGS}"
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
        ARCHIVE_OUTPUT_DIRECTORY "${READDY_PLUGIN_OUTPUT_DIRECTORY}"
        LIBRARY_OUTPUT_DIRECTORY "${READDY_PLUGIN_OUTPUT_DIRECTORY}"
        LINKER_LANGUAGE CXX)

INSTALL(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION readdy/readdy_plugins
        LIBRARY DESTINATION readdy/readdy_plugins
        ARCHIVE DESTINATION readdy/readdy_plugins
        )

IF (READDY_CREATE_TEST_TARGET)
    ADD_SUBDIRECTORY(test)
ENDIF()