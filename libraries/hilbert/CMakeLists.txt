project(hilbert C)

set(HILBERT_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" CACHE INTERNAL "Hilbert include directory")
set(HILBERT_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/hilbert.c")

add_library(${PROJECT_NAME} SHARED ${HILBERT_SRC})
target_include_directories(${PROJECT_NAME} PUBLIC ${HILBERT_INCLUDE_DIR})

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
