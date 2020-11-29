set(SOURCES
    QASM_to_Sylvan.c
)

add_executable(QASM_to_Sylvan QASM_to_Sylvan.c ${SOURCES})
target_link_libraries(QASM_to_Sylvan sylvan)

# Include profiler if installed
include(CheckIncludeFiles)
check_include_files("gperftools/profiler.h" HAVE_PROFILER)

if(HAVE_PROFILER)
    set_target_properties(bench_algs PROPERTIES COMPILE_DEFINITIONS "HAVE_PROFILER")
    target_link_libraries(bench_algs profiler)
endif()
