SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/cpu/src")
SET(CPU_INCLUDE_DIR "${READDY_GLOBAL_DIR}/kernels/cpu/include")

# --- main sources ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/CPUKernel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/model/CPUNeighborList.cpp")

# --- programs ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/CPUProgramFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/CPUEulerBDIntegrator.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/CPUDefaultReactionProgram.cpp")
