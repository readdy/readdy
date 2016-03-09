#INCLUDE(ExternalProject)
SET(Boost.DLL_ROOT "${READDY_GLOBAL_DIR}/libraries/dll")
SET(Boost.DLL_INCLUDE_DIRS "${Boost.DLL_ROOT}/include")
#ExternalProject_Add(Boost.DLL
#        DOWNLOAD_COMMAND ""
#        UPDATE_COMMAND ""
#        PATCH_COMMAND ""
#
#        SOURCE_DIR ${Boost.DLL_ROOT}
#        INSTALL_DIR ${Boost.DLL_ROOT}
#        BUILD_IN_SOURCE 1
#        INSTALL_COMMAND ""
#        )