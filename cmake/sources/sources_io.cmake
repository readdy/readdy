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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/readdy/main/io")

SET(H5RD_INCLUDE_DIR "${READDY_GLOBAL_DIR}/libraries/h5rd/include")

# hdf5
# SET(HDF5_USE_STATIC_LIBRARIES OFF)
FIND_PACKAGE(HDF5 COMPONENTS HL REQUIRED)
MESSAGE(STATUS "HDF5_FOUND: ${HDF5_FOUND}, version ${HDF5_VERSION}")

# includes
SET(BLOSC_INCLUDE_DIR "${READDY_GLOBAL_DIR}/libraries/c-blosc/blosc")
MESSAGE(STATUS "BLOSC INCLUDE DIR: ${BLOSC_INCLUDE_DIR}")
SET(IO_INCLUDE_DIRS "${COMMON_INCLUDE_DIRS};${HDF5_INCLUDE_DIRS};${HDF5_HL_INCLUDE_DIR};${BLOSC_INCLUDE_DIR};${H5RD_INCLUDE_DIR}" CACHE INTERNAL "IO Include dirs")

# libraries
SET(BLOSC_LIBRARIES "blosc_shared")
SET(READDY_IO_LIBRARIES "${READDY_COMMON_LIBRARIES};${BLOSC_LIBRARIES};${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}" CACHE INTERNAL "IO Libraries")
LIST(REMOVE_DUPLICATES READDY_IO_LIBRARIES)

# sources
LIST(APPEND READDY_IO_SOURCES "${SOURCES_DIR}/BloscFilter.cpp")
LIST(APPEND READDY_IO_SOURCES "${SOURCES_DIR}/blosc_filter.h")
LIST(APPEND READDY_IO_SOURCES "${SOURCES_DIR}/blosc_filter.c")
LIST(APPEND READDY_IO_SOURCES "${SOURCES_DIR}/blosc_filter.h")
LIST(APPEND READDY_IO_SOURCES "${SOURCES_DIR}/blosc_plugin.c")

# all sources
LIST(APPEND READDY_ALL_SOURCES ${READDY_IO_SOURCES})