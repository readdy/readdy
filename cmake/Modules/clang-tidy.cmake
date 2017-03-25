find_program(CLANG_TIDY "clang-tidy")
if(CLANG_TIDY)
    macro(add_clang_tidy_target TARGET_NAME SOURCE_FILES INCLUDES)
        message(status "added clang-tidy as ${TARGET_NAME}" with source_files=[${SOURCE_FILES}] and includes=[${INCLUDES}])
        foreach(inc ${INCLUDES})
            list(APPEND ALL_INCLUDES "-I${inc}")
        endforeach()
        message(status "got includes ${ALL_INCLUDES}")
        add_custom_target(${TARGET_NAME}
                COMMAND ${CLANG_TIDY}
                -p ${CMAKE_SOURCE_DIR}
                --header-filter=${CMAKE_SOURCE_DIR}/include
                --checks="*"
                ${SOURCE_FILES}
                --
                ${ALL_INCLUDES})
    endmacro()
else()
    message(STATUS "Could not find clang-tidy.")
endif()