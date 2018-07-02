cmake_minimum_required(VERSION 3.11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

find_package(READDY REQUIRED)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../libraries/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/pybind11)

pybind11_add_module(custom_integrator_example binding.cpp)
target_include_directories(custom_integrator_example PUBLIC "${READDY_INCLUDE_DIRS}")
target_link_libraries(custom_integrator_example LINK_PUBLIC ${READDY_LIBRARIES})
