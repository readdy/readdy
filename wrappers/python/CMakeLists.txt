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


# import python
MESSAGE(STATUS "Python library set: ${PYTHON_LIBRARY}")

# set output directory variable for modules
IF(READDY_DEBUG_PYTHON_MODULES)
    SET(PYTHON_CPP_MODULES_OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/python/readdy/_internal")
ELSE()
    SET(PYTHON_CPP_MODULES_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/readdy/_internal")
ENDIF()

# add modules
ADD_SUBDIRECTORY(src/cxx)

# install via setup.py
FIND_PROGRAM(PYTHON "python")
if (PYTHON)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/src/python/setup.py")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(VERSION_IN  "${CMAKE_CURRENT_SOURCE_DIR}/src/python/readdy/_version.py.in")
    set(VERSION_OUT "${CMAKE_CURRENT_SOURCE_DIR}/src/python/readdy/_version.py")
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})
    configure_file(${VERSION_IN} ${VERSION_OUT})

    add_custom_command(OUTPUT ${OUTPUT}
            COMMAND ${PYTHON} ${SETUP_PY} build
            COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
            DEPENDS readdybinding)

    add_custom_target(target ALL DEPENDS ${OUTPUT})

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX})")
ELSE()
    MESSAGE(FATAL_ERROR "No python executable was found!")
endif()