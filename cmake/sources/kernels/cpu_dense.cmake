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


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/cpu_dense/src")
SET(CPU_DENSE_INCLUDE_DIR "${READDY_GLOBAL_DIR}/kernels/cpu_dense/include")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/CPUDKernel.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/CPUDStateModel.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/model/CPUDParticleData.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/model/CPUDNeighborList.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/observables/CPUDObservableFactory.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/observables/CPUDObservables.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/potentials/CPUDPotentialFactory.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/actions/reactions/Event.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/actions/reactions/CPUDGillespie.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/actions/reactions/CPUDGillespieParallel.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/actions/reactions/ReactionUtils.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/actions/reactions/CPUDUncontrolledApproximation.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/actions/CPUDCompartments.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/actions/CPUDEulerBDIntegrator.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/actions/CPUDActionFactory.cpp")
