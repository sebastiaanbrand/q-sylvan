set(SOURCES
    grover.c
    supremacy.c
)

add_executable(bench bench.c ${SOURCES})
target_link_libraries(bench sylvan)

add_executable(test_algs test_algs.c ${SOURCES})
target_link_libraries(test_algs sylvan)