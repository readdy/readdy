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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/readdy/main/model")

# libraries
SET(READDY_MODEL_LIBRARIES "${READDY_IO_LIBRARIES}" CACHE INTERNAL "Model libraries")

# includes
SET(MODEL_INCLUDE_DIRS "${IO_INCLUDE_DIRS}" CACHE INTERNAL "Model include dirs")

# sources
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/KernelContext.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Particle.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Vec3.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Kernel.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Observables.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Aggregators.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Potentials.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Trajectory.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Actions.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Reactions.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Utils.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/IOUtils.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/Compartments.cpp")

# topologies
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/Topology.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/TopologyPotential.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/BondedPotential.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/AnglePotential.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/topologies/TorsionPotential.cpp")

# observables
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/CenterOfMass.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/HistogramAlongAxis.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Particles.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/NParticles.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Forces.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Reactions.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/ReactionCounts.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/Positions.cpp")
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/observables/RadialDistribution.cpp")

# internal sources
LIST(APPEND READDY_MODEL_SOURCES "${SOURCES_DIR}/_internal/ObservableWrapper.cpp")

# all sources
LIST(APPEND READDY_ALL_SOURCES ${READDY_MODEL_SOURCES})