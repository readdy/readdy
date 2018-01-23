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