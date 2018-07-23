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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/readdy/main/model")

# libraries
SET(READDY_MODEL_LIBRARIES "${READDY_IO_LIBRARIES}" CACHE INTERNAL "Model libraries")

# includes
SET(MODEL_INCLUDE_DIRS "${IO_INCLUDE_DIRS}" CACHE INTERNAL "Model include dirs")

# sources
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Context.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/StateModel.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Particle.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/ParticleTypeRegistry.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Aggregators.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Actions.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Utils.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/IOUtils.cpp")

# compartments
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/compartments/Compartments.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/compartments/CompartmentRegistry.cpp")

# potentials
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/potentials/Potentials.cpp")

# reactions
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/reactions/Reactions.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/reactions/ReactionRegistry.cpp")

# topologies
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/Utils.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/GraphTopology.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/TopologyRegistry.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/graph/Graph.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/potentials/BondedPotential.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/potentials/AnglePotential.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/potentials/TorsionPotential.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/reactions/Recipe.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/reactions/TopologyReactionActions.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/reactions/StructuralTopologyReaction.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/reactions/SpatialTopologyReaction.cpp")

# observables
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Types.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Trajectory.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/HistogramAlongAxis.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Particles.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/NParticles.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Forces.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Reactions.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/ReactionCounts.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Energy.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Positions.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Topologies.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/RadialDistribution.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Virial.cpp")

# all sources
LIST(APPEND READDY_ALL_SOURCES ${READDY_MODEL_SOURCES})