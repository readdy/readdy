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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/cpu/src")
SET(CPU_INCLUDE_DIR "${READDY_GLOBAL_DIR}/kernels/cpu/include")

# --- main sources ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/CPUKernel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/CPUStateModel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/observables/CPUObservableFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/observables/CPUObservables.cpp")

# --- neighbor list ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/nl/CellLinkedList.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/nl/ContiguousCellLinkedList.cpp")

# --- actions ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/CPUActionFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/CPUEulerBDIntegrator.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/CPUCalculateForces.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/CPUEvaluateCompartments.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/CPUEvaluateTopologyReactions.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/reactions/ReactionUtils.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/reactions/Event.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/reactions/CPUUncontrolledApproximation.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/reactions/CPUGillespie.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/topologies/CPUTopologyActions.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/actions/topologies/CPUTopologyActionFactory.cpp")

LIST(APPEND READDY_ALL_SOURCES ${CPU_SOURCES})
