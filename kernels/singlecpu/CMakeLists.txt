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


PROJECT(readdy_kernel_singlecpu CXX)

# sources and headers
INCLUDE("${READDY_GLOBAL_DIR}/cmake/sources/kernels/singlecpu.cmake")

# create library
ADD_LIBRARY(${PROJECT_NAME} SHARED ${SINGLECPU_SOURCES} ${READDY_INCLUDE_DIRS})

# includes
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS})

# link
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CMAKE_DL_LIBS} ${READDY_SINGLECPU_DEPENDENT_LIBRARIES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} readdy_common readdy_model)

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    LINK_FLAGS "${EXTRA_LINK_FLAGS}"
    COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
)

INSTALL(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)