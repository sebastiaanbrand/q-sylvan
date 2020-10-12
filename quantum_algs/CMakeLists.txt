set(SOURCES
    grover.c
    random_circuit.c
    shor.c
    supremacy.c
)

add_executable(bench_algs bench_algs.c ${SOURCES})
target_link_libraries(bench_algs sylvan)

add_executable(test_algs test_algs.c ${SOURCES})
target_link_libraries(test_algs sylvan)

# Include profiler if installed
include(CheckIncludeFiles)
check_include_files("gperftools/profiler.h" HAVE_PROFILER)

if(HAVE_PROFILER)
    set_target_properties(bench_algs PROPERTIES COMPILE_DEFINITIONS "HAVE_PROFILER")
    target_link_libraries(bench_algs profiler)
endif()
