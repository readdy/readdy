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