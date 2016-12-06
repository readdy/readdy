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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/singlecpu/src")

# libraries
SET(READDY_SINGLECPU_DEPENDENT_LIBRARIES "${READDY_COMMON_LIBRARIES};${READDY_MODEL_LIBRARIES};${READDY_IO_LIBRARIES}" CACHE INTERNAL "Readdy libraries")
LIST(REMOVE_DUPLICATES READDY_SINGLECPU_DEPENDENT_LIBRARIES)

# sources
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/SCPUKernel.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/SCPUProgramFactory.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/SCPUStateModel.cpp")

# --- programs ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SCPUTestProgram.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SCPUAddParticle.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SCPUEulerBDIntegrator.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SCPUUpdateNeighborList.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SCPUCalculateForces.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SCPUReactionImpls.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SCPUCompartments.cpp")

# --- potentials ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/potentials/SCPUPotentialFactory.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/potentials/SCPUPotentialsOrder1.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/potentials/SCPUPotentialsOrder2.cpp")

# --- reactions ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/reactions/SCPUReactionFactory.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/reactions/SCPUReactions.cpp")

# --- model ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/model/SCPUParticleData.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/model/SCPUNeighborList.cpp")

# --- observables ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/observables/SCPUObservableFactory.cpp")


# --- all sources ---
LIST(APPEND READDY_ALL_SOURCES ${SINGLECPU_SOURCES})