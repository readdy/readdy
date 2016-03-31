#INCLUDE(ExternalProject)
SET(Boost.DLL_ROOT "${READDY_GLOBAL_DIR}/libraries/dll")
SET(Boost.DLL_INCLUDE_DIRS "${Boost.DLL_ROOT}/include")

INSTALL(DIRECTORY "${Boost.DLL_INCLUDE_DIRS}/boost" COMPONENT Development DESTINATION include)