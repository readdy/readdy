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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/cpu/src")
SET(CPU_INCLUDE_DIR "${READDY_GLOBAL_DIR}/kernels/cpu/include")

# --- main sources ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/CPUKernel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/CPUStateModel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/model/CPUNeighborList.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/model/CPUParticleData.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/potentials/CPUPotentialFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/observables/CPUObservableFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/observables/CPUObservables.cpp")

# --- programs ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/CPUProgramFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/CPUEulerBDIntegrator.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/CPUCompartments.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/ReactionUtils.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/Event.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/CPUUncontrolledApproximation.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/CPUGillespie.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/CPUGillespieParallel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/FilteredGillespieParallel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/NextSubvolumesReactionScheduler.cpp")

# --- c sources ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/util/hilbert.c")
