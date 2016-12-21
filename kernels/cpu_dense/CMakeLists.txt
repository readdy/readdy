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


PROJECT(ReaDDy_kernel_cpu_dense C CXX)

# sources and headers
INCLUDE("${READDY_GLOBAL_DIR}/cmake/sources/kernels/cpu_dense.cmake")

# create library
ADD_LIBRARY(${PROJECT_NAME} SHARED ${CPU_DENSE_SOURCES} ${ReaDDy_INCLUDE_DIRS} ${CPU_DENSE_INCLUDE_DIR})

# includes
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC ${ReaDDy_INCLUDE_DIRS} ${CPU_DENSE_INCLUDE_DIR})

# link
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CMAKE_DL_LIBS} ReaDDy)
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
        LINK_FLAGS "${EXTRA_LINK_FLAGS}"
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
        CXX_STANDARD ${READDY_STANDARD}
        CXX_STANDARD_REQUIRED ON
        )

INSTALL(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION readdy/readdy_plugins
        LIBRARY DESTINATION readdy/readdy_plugins
        ARCHIVE DESTINATION readdy/readdy_plugins
        )

IF (READDY_CREATE_TEST_TARGET)
    ADD_SUBDIRECTORY(test)
ENDIF()