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

find_path(READDY_INCLUDE_DIR NAMES readdy/readdy.h DOC "The ReaDDy include directory")
find_library(READDY_LIBRARY NAMES readdy DOC "The ReaDDy library")
find_library(READDY_CPU_LIBRARY NAMES readdy_kernel_cpu
        HINTS "${CMAKE_PREFIX_PATH}/readdy/readdy_plugins"
        DOC "The ReaDDy CPU kernel")
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(READDY REQUIRED_VARS READDY_LIBRARY READDY_INCLUDE_DIR READDY_CPU_LIBRARY)
if(READDY_FOUND)
    set(READDY_LIBRARIES ${READDY_LIBRARY} ${READDY_CPU_LIBRARY})
    set(READDY_INCLUDE_DIRS ${READDY_INCLUDE_DIR})
    if(NOT TARGET READDY::READDY)
        add_library(READDY::READDY UNKNOWN IMPORTED)
        set_target_properties(READDY::READDY PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${READDY_INCLUDE_DIR}")
        set_property(TARGET READDY::READDY APPEND PROPERTY IMPORTED_LOCATION "${READDY_LIBRARY}")
    endif()
    if(NOT TARGET READDY::READDY_CPU)
        add_library(READDY::READDY_CPU UNKNOWN IMPORTED)
        set_target_properties(READDY::READDY_CPU PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${READDY_INCLUDE_DIR}")
        set_property(TARGET READDY::READDY_CPU APPEND PROPERTY IMPORTED_LOCATION "${READDY_CPU_LIBRARY}")
    endif()
endif()

mark_as_advanced(READDY_INCLUDE_DIR READDY_LIBRARY)
