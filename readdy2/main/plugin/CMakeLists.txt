PROJECT(ReaDDy2_plugin)

# requires boost dll
INCLUDE(ExternalProject)
SET(Boost.DLL_ROOT "${READDY_GLOBAL_DIR}/libraries/external/dll")
SET(Boost.DLL_INCLUDE_DIRS "${Boost.DLL_ROOT}/include/boost")
ExternalProject_Add(Boost.DLL
        PREFIX ${Boost.DLL_ROOT}
        URL https://github.com/boostorg/dll/archive/boost-review.tar.gz
        BINARY_DIR ${Boost.DLL_ROOT}
        INSTALL_DIR ${Boost.DLL_ROOT}
        INSTALL_COURLMMAND ""
)

# sources
LIST(APPEND READDY_PLUGIN_SOURCES
        KernelProvider.cpp)

# expose includes
SET(${PROJECT_NAME}_INCLUDE_DIRS "${ReaDDy2_common_INCLUDE_DIRS};${Boost.DLL_INCLUDE_DIRS};${READDY_GLOBAL_INCLUDE_DIR}/plugin" CACHE INTERNAL "${PROJECT_NAME}: Include directories" FORCE)

# library
ADD_LIBRARY(${PROJECT_NAME} ${READDY_PLUGIN_SOURCES} "${${PROJECT_NAME}_INCLUDE_DIRS}")

# includes
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC ${${PROJECT_NAME}_INCLUDE_DIRS})

# Required on Unix OS family to be able to be linked into shared libraries.
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# link
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ReaDDy2_common)

# use c++11
set_property (TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS -std=c++11)