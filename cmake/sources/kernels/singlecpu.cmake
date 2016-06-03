SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/singlecpu")

# sources
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/SingleCPUKernel.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/SingleCPUProgramFactory.cpp")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/SingleCPUKernelStateModel.cpp")

# headers
LIST(APPEND SINGLECPU_HEADERS "${SOURCES_DIR}/SingleCPUKernel.h")
LIST(APPEND SINGLECPU_HEADERS "${SOURCES_DIR}/SingleCPUProgramFactory.h")
LIST(APPEND SINGLECPU_HEADERS "${SOURCES_DIR}/SingleCPUKernelStateModel.h")

# --- programs ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SingleCPUTestProgram.cpp")
LIST(APPEND SINGLECPU_HEADERS "${SOURCES_DIR}/programs/SingleCPUTestProgram.h")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SingleCPUAddParticleProgram.cpp")
LIST(APPEND SINGLECPU_HEADERS "${SOURCES_DIR}/programs/SingleCPUAddParticleProgram.h")
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/programs/SingleCPUDiffuseProgram.cpp")
LIST(APPEND SINGLECPU_HEADERS "${SOURCES_DIR}/programs/SingleCPUDiffuseProgram.h")

# --- potentials ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/potentials/PotentialsOrder1.cpp")
LIST(APPEND SINGLECPU_HEADERS "${SOURCES_DIR}/potentials/PotentialsOrder1.h")

# --- model ---
LIST(APPEND SINGLECPU_SOURCES "${SOURCES_DIR}/model/SingleCPUParticleData.cpp")
LIST(APPEND SINGLECPU_HEADERS "${SOURCES_DIR}/model/SingleCPUParticleData.h")

