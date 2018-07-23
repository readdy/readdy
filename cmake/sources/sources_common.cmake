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
