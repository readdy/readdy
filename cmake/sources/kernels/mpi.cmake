####################################################################
# Copyright © 2019 Computational Molecular Biology Group,          #
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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/mpi/src")
SET(MPI_INCLUDE_DIR "${READDY_GLOBAL_DIR}/kernels/mpi/include")

# sources
LIST(APPEND MPI_SOURCES "${SOURCES_DIR}/MPIKernel.cpp")
LIST(APPEND MPI_SOURCES "${SOURCES_DIR}/MPIStateModel.cpp")
LIST(APPEND MPI_SOURCES "${SOURCES_DIR}/actions/MPIActionFactory.cpp")

# --- actions ---
LIST(APPEND MPI_SOURCES "${SOURCES_DIR}/actions/MPICreateNeighborList.cpp")
LIST(APPEND MPI_SOURCES "${SOURCES_DIR}/actions/MPIReactionImpls.cpp")
#LIST(APPEND MPI_SOURCES "${SOURCES_DIR}/actions/MPIEvaluateCompartments.cpp")
#LIST(APPEND MPI_SOURCES "${SOURCES_DIR}/actions/MPIEvaluateTopologyReactions.cpp")

# --- model ---
#LIST(APPEND MPI_SOURCES "${SOURCES_DIR}/MPIParticleData.cpp")

# --- observables ---
LIST(APPEND MPI_SOURCES "${SOURCES_DIR}/MPIObservableFactory.cpp")

# --- topology actions ---
#LIST(APPEND MPI_SOURCES "${SOURCES_DIR}/topologies/MPITopologyActionFactory.cpp")

# --- all sources ---
LIST(APPEND READDY_ALL_SOURCES ${MPI_SOURCES})
