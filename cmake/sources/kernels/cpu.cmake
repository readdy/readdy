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
