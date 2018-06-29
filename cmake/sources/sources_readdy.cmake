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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/readdy/main")

# includes
SET(READDY_INCLUDE_DIRS "${READDY_GLOBAL_INCLUDE_DIR};${IO_INCLUDE_DIRS};${PLUGIN_INCLUDE_DIRS}" CACHE INTERNAL "Readdy include dirs")
LIST(REMOVE_DUPLICATES READDY_INCLUDE_DIRS)

# libraries
SET(READDY_DEPENDENT_LIBRARIES "${READDY_COMMON_LIBRARIES};${READDY_MODEL_LIBRARIES};${READDY_IO_LIBRARIES};${READDY_PLUGIN_LIBRARIES};${READDY_SINGLECPU_DEPENDENT_LIBRARIES}" CACHE INTERNAL "Readdy libraries")
LIST(REMOVE_DUPLICATES READDY_DEPENDENT_LIBRARIES)

# sources
LIST(APPEND READDY_MAIN_SOURCES "")

# all sources
LIST(APPEND READDY_ALL_SOURCES ${READDY_MAIN_SOURCES})