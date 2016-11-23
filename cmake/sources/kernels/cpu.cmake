SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/cpu/src")
SET(CPU_INCLUDE_DIR "${READDY_GLOBAL_DIR}/kernels/cpu/include")

# --- main sources ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/Kernel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/StateModel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/model/NeighborList.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/model/ParticleData.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/potentials/PotentialFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/observables/ObservableFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/observables/Observables.cpp")

# --- programs ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/ProgramFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/EulerBDIntegrator.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/Compartments.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/ReactionUtils.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/Event.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/UncontrolledApproximation.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/Gillespie.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/GillespieParallel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/FilteredGillespieParallel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/reactions/NextSubvolumesReactionScheduler.cpp")

# --- c sources ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/util/hilbert.c")
