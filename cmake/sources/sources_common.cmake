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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/readdy/main/common")

SET(SPDLOG_INCLUDE_DIR "${READDY_GLOBAL_DIR}/libraries/spdlog/include")
SET(JSON_INCLUDE_DIR "${READDY_GLOBAL_DIR}/libraries/json/include")

FIND_PACKAGE(Threads REQUIRED)

# include dirs
SET(COMMON_INCLUDE_DIRS "${READDY_GLOBAL_INCLUDE_DIR};${SPDLOG_INCLUDE_DIR};${JSON_INCLUDE_DIR}" CACHE INTERNAL "Common include dirs")

# dependent libraries
SET(READDY_COMMON_LIBRARIES "" CACHE INTERNAL "Common libraries")

# sources
LIST(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/../api/KernelConfiguration.cpp")
LIST(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/Utils.cpp")
LIST(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/filesystem.cpp")
LIST(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/Config.cpp")
LIST(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/logging.cpp")
LIST(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/Timer.cpp")

# all sources
LIST(APPEND READDY_ALL_SOURCES ${READDY_COMMON_SOURCES})
