set(SOURCES
    QASM_to_circuit.c
    circuit.c
)

add_executable(test_circuit test_circuit.c ${SOURCES})
target_link_libraries(test_circuit qsylvan)
add_executable(circuit_to_Sylvan circuit_to_Sylvan.c ${SOURCES})
target_link_libraries(circuit_to_Sylvan qsylvan popt)

