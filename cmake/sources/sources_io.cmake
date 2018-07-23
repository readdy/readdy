####################################################################
# Copyright © 2018 Computational Molecular Biology Group,          #
#                  Freie Universität Berlin (GER)                  #
#                                                                  #
# Redistribution and use in source and binary forms, with or       #
# without modification, are permitted provided that the            #
# following conditions are met:                                    #
#  1. Redistributions of source code must retain the above         #
#     copyright notice, this list of conditions and the            #
#     following disclaimer.                                        #
#  2. Redistributions in binary form must reproduce the above      #
#     copyright notice, this list of conditions and the following  #
#     disclaimer in the documentation and/or other materials       #
#     provided with the distribution.                              #
#  3. Neither the name of the copyright holder nor the names of    #
#     its contributors may be used to endorse or promote products  #
#     derived from this software without specific                  #
#     prior written permission.                                    #
#                                                                  #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           #
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      #
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         #
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         #
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            #
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     #
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         #
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER #
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      #
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       #
####################################################################


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