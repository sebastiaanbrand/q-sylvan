set(SOURCES
    grover_cnf.c
    grover.c
    random_circuit.c
    shor.c
    supremacy.c
)

# Bell state example
add_executable(bell_state bell_state.c)
target_link_libraries(bell_state qsylvan)

# Executable to run algorithms
add_executable(alg_run alg_run.c ${SOURCES})
target_link_libraries(alg_run qsylvan) # argp

# Tests
add_executable(test_algs test_algs.c ${SOURCES})
target_link_libraries(test_algs qsylvan)
add_test(test_algs test_algs)

