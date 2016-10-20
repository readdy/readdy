SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/singlecpu/src")

# libraries
SET(READDY_SINGLECPU_DEPENDENT_LIBRARIES "${READDY_COMMON_LIBRARIES};${READDY_MODEL_LIBRARIES};${READDY_IO_LIBRARIES}" CACHE INTERNAL "Readdy libraries")
LIST(REMOVE_DUPLICATES READDY_SINGLECPU_DEPENDENT_LIBRARIES)

# sources
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/SingleCPUKernel.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/SingleCPUProgramFactory.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/SingleCPUKernelStateModel.cpp")

# --- programs ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SingleCPUTestProgram.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SingleCPUAddParticleProgram.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SingleCPUEulerBDIntegrator.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SingleCPUUpdateNeighborList.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SingleCPUCalculateForces.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SingleCPUReactionImpls.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/Compartments.cpp")

# --- potentials ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/potentials/SingleCPUPotentialFactory.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/potentials/PotentialsOrder1.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/potentials/PotentialsOrder2.cpp")

# --- reactions ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/reactions/SingleCPUReactionFactory.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/reactions/SingleCPUReactions.cpp")

# --- model ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/model/SingleCPUParticleData.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/model/SingleCPUNeighborList.cpp")

# --- observables ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/observables/SingleCPUObservableFactory.cpp")


# --- all sources ---
LIST(APPEND READDY_ALL_SOURCES ${SINGLECPU_SOURCES})