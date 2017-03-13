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


LIST(APPEND API_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/PySimulation.cpp)
LIST(APPEND API_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/ExportSchemeApi.h)
LIST(APPEND API_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/ExportTopologies.cpp)
LIST(APPEND API_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/ExportObservables.h)
LIST(APPEND API_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/PyFunction.h)
LIST(APPEND API_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/PyPotential.h)
LIST(APPEND API_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/PyPotential.cpp)

PYBIND11_ADD_MODULE(api ${API_SOURCES})
TARGET_INCLUDE_DIRECTORIES(api PUBLIC "${READDY_INCLUDE_DIRS}")
TARGET_LINK_LIBRARIES(api LINK_PUBLIC ${PYTHON_LIBRARIES} ${READDY_DEPENDENT_LIBRARIES} readdy)

SET_TARGET_PROPERTIES(api
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${PYTHON_CPP_MODULES_OUTPUT_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${PYTHON_CPP_MODULES_OUTPUT_DIR}")