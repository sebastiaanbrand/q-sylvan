set(SOURCES
    QASM_to_circuit.c
    circuit.c
)

add_executable(test_circuit test_circuit.c ${SOURCES})
target_link_libraries(test_circuit qsylvan)
add_executable(qsylvan_qasm qsylvan_qasm.c ${SOURCES})
target_link_libraries(qsylvan_qasm qsylvan popt)

add_executable(sim_qasm sim_qasm.c simple_parser.cpp)
target_link_libraries(sim_qasm qsylvan)

