project(readdy_examples)

add_executable(diffusion_benchmark DiffusionBenchmark.cpp)
target_include_directories(diffusion_benchmark PUBLIC ${READDY_INCLUDE_DIRS})
target_link_libraries(diffusion_benchmark PUBLIC readdy)
set_target_properties(diffusion_benchmark PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
install(TARGETS diffusion_benchmark
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)