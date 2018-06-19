cmake_minimum_required(VERSION 3.11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

find_package(READDY REQUIRED)

message(STATUS "current source dir ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "source dir ${CMAKE_SOURCE_DIR}")
if("${CMAKE_CURRENT_SOURCE_DIR}" EQUAL "${CMAKE_SOURCE_DIR}")
    message(STATUS "adding pybind subdir")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../libraries/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/pybind11)
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../libraries/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/pybind11)

set(SRC binding.cpp TestPotential.h TestPotential.cpp)
pybind11_add_module(custom_potential_example binding.cpp)
target_include_directories(custom_potential_example PUBLIC "${READDY_INCLUDE_DIRS}")
target_link_libraries(custom_potential_example LINK_PUBLIC ${READDY_LIBRARIES})
