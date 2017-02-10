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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/readdy/main/plugin")

# includes
SET(PLUGIN_INCLUDE_DIRS "${MODEL_INCLUDE_DIRS}" CACHE INTERNAL "Plugin include dirs")

# libraries
SET(READDY_PLUGIN_LIBRARIES "${READDY_MODEL_LIBRARIES};${CMAKE_DL_LIBS};${READDY_SINGLECPU_DEPENDENT_LIBRARIES}" CACHE INTERNAL "Plugin libraries")

# sources
LIST(APPEND READDY_PLUGIN_SOURCES "${SOURCES_DIR}/KernelProvider.cpp")

# all sources
LIST(APPEND READDY_ALL_SOURCES ${READDY_PLUGIN_SOURCES})