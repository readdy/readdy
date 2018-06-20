#####################################################################
# Copyright (c) 2018 Computational Molecular Biology Group,         #
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

find_path(READDY_INCLUDE_DIR NAMES readdy/readdy.h DOC "The ReaDDy include directory")
find_library(READDY_LIBRARY NAMES readdy DOC "The ReaDDy library")
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(READDY REQUIRED_VARS READDY_LIBRARY READDY_INCLUDE_DIR)
if(READDY_FOUND)
    set( READDY_LIBRARIES ${READDY_LIBRARY} )
    set( READDY_INCLUDE_DIRS ${READDY_INCLUDE_DIR} )
    if(NOT TARGET READDY::READDY)
        add_library(READDY::READDY UNKNOWN IMPORTED)
        set_target_properties(READDY::READDY PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${READDY_INCLUDE_DIR}")
        set_property(TARGET READDY::READDY APPEND PROPERTY IMPORTED_LOCATION "${READDY_LIBRARY}")
    endif()
endif()

mark_as_advanced(READDY_INCLUDE_DIR READDY_LIBRARY)
