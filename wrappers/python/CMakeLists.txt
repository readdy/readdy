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

cmake_minimum_required(VERSION 3.11)

project(readdy_python_wrapper C CXX)

# import python
MESSAGE(STATUS "Python library set: ${PYTHON_LIBRARY}")

# set output directory variable for modules
if(READDY_DEBUG_PYTHON_MODULES)
    SET(PYTHON_CPP_MODULES_OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/python/readdy/_internal")
else()
    SET(PYTHON_CPP_MODULES_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/readdy/_internal")
endif()

add_custom_target(CopyPythonSources
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/src/python/readdy ${CMAKE_CURRENT_BINARY_DIR}/readdy)
# add modules
ADD_SUBDIRECTORY(src/cxx)

install (
        DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/src/python/readdy/
        DESTINATION ${INSTALL_DIR}
        FILES_MATCHING PATTERN "*.py"
)
configure_file(${CMAKE_CURRENT_LIST_DIR}/src/python/readdy/_version.py.in ${CMAKE_CURRENT_BINARY_DIR}/src/python/readdy/_version.py)

install (FILES ${CMAKE_CURRENT_BINARY_DIR}/src/python/readdy/_version.py
         DESTINATION ${INSTALL_DIR}
)
