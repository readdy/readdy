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
