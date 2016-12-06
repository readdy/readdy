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


LIST(APPEND COMMON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CommonModule.cpp)

PYBIND11_ADD_MODULE(common ${COMMON_SOURCES})
TARGET_INCLUDE_DIRECTORIES(common PUBLIC "${ReaDDy_INCLUDE_DIRS}")
TARGET_LINK_LIBRARIES(common LINK_PUBLIC ${READDY_DEPENDENT_LIBRARIES} ReaDDy)

SET_TARGET_PROPERTIES(common
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${PYTHON_CPP_MODULES_OUTPUT_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${PYTHON_CPP_MODULES_OUTPUT_DIR}"
        CXX_STANDARD ${READDY_STANDARD}
        CXX_STANDARD_REQUIRED ON
        )