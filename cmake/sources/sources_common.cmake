set(SOURCES_DIR "${READDY_GLOBAL_DIR}/readdy/main/common")

find_package(Threads REQUIRED)

# include dirs
set(COMMON_INCLUDE_DIRS "${READDY_GLOBAL_INCLUDE_DIR};${blosc_INCLUDE_DIR}" CACHE INTERNAL "Common include dirs")

# dependent libraries
set(READDY_COMMON_LIBRARIES "nlohmann_json::nlohmann_json;spdlog::spdlog_header_only;blosc::blosc;fmt::fmt-header-only" CACHE INTERNAL "Common libraries")

# sources
list(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/../api/KernelConfiguration.cpp")
list(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/Utils.cpp")
list(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/filesystem.cpp")
list(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/Config.cpp")
list(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/logging.cpp")
list(APPEND READDY_COMMON_SOURCES "${SOURCES_DIR}/Timer.cpp")

# all sources
list(APPEND READDY_ALL_SOURCES ${READDY_COMMON_SOURCES})
