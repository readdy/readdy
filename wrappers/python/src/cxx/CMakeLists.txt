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


# add READDY_WITH_PYTHON define
ADD_DEFINITIONS(-DREADDY_WITH_PYTHON)

set(SOURCE_FILES
        ReaddyBinding.cpp
        common/common.h
        common/CommonModule.cpp
        common/ExportIO.cpp
        common/ReadableReactionRecord.h
        common/ReadableReactionRecord.cpp
        common/ReadableParticle.h
        common/SpdlogPythonSink.h
        common/Utils.cpp
        api/ExportObservables.h
        api/ExportSchemeApi.cpp
        api/ExportTopologies.cpp
        api/ExportKernelContext.cpp
        api/PyFunction.h
        api/ApiModule.cpp)
pybind11_add_module(readdybinding ${SOURCE_FILES})
target_include_directories(readdybinding PUBLIC "${READDY_INCLUDE_DIRS}")
target_link_libraries(readdybinding LINK_PUBLIC ${READDY_DEPENDENT_LIBRARIES} readdy)

SET_TARGET_PROPERTIES(readdybinding
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${PYTHON_CPP_MODULES_OUTPUT_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${PYTHON_CPP_MODULES_OUTPUT_DIR}")
