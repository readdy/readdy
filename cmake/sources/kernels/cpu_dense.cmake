SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/cpu_dense/src")
SET(CPU_DENSE_INCLUDE_DIR "${READDY_GLOBAL_DIR}/kernels/cpu_dense/include")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/CPUDKernel.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/CPUDStateModel.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/model/CPUDParticleData.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/model/CPUDNeighborList.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/observables/CPUDObservableFactory.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/observables/CPUDObservables.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/potentials/CPUDPotentialFactory.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/reactions/Event.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/reactions/CPUDGillespie.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/reactions/CPUDGillespieParallel.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/reactions/ReactionUtils.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/reactions/CPUDUncontrolledApproximation.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/CPUDCompartments.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/CPUDEulerBDIntegrator.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/CPUDProgramFactory.cpp")
