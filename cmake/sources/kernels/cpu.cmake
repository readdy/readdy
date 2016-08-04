SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/cpu/src")
SET(CPU_INCLUDE_DIR "${READDY_GLOBAL_DIR}/kernels/cpu/include")

# --- main sources ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/CPUKernel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/CPUStateModel.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/potentials/CPUPotentialFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/observables/ObservableFactory.cpp")

# --- programs ---
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/CPUProgramFactory.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/CPUEulerBDIntegrator.cpp")
LIST(APPEND CPU_SOURCES "${SOURCES_DIR}/programs/Reactions.cpp")
