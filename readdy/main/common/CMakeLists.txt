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


# general utilities.
PROJECT(readdy_common CXX)

INCLUDE("${READDY_GLOBAL_DIR}/cmake/sources/sources_common.cmake")

ADD_LIBRARY(${PROJECT_NAME} ${READDY_COMMON_SOURCES} "${READDY_GLOBAL_INCLUDE_DIR}/readdy/common")

TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC ${COMMON_INCLUDE_DIRS})

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${EXTRA_LINK_FLAGS}"
    COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${READDY_COMMON_LIBRARIES})

INSTALL(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)