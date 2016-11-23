SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/cpu_dense/src")
SET(CPU_DENSE_INCLUDE_DIR "${READDY_GLOBAL_DIR}/kernels/cpu_dense/include")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/Kernel.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/StateModel.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/model/ParticleData.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/model/NeighborList.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/observables/ObservableFactory.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/observables/Observables.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/potentials/PotentialFactory.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/reactions/Event.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/reactions/Gillespie.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/reactions/GillespieParallel.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/reactions/ReactionUtils.cpp")

LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/Compartments.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/EulerBDIntegrator.cpp")
LIST(APPEND CPU_DENSE_SOURCES "${SOURCES_DIR}/programs/ProgramFactory.cpp")
