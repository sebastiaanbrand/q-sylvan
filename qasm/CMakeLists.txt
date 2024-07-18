#set(SOURCES
#    QASM_to_circuit.c
#    circuit.c
#)
#add_executable(test_circuit test_circuit.c ${SOURCES})
#target_link_libraries(test_circuit qsylvan)
#add_executable(qsylvan_qasm qsylvan_qasm.c ${SOURCES})
#target_link_libraries(qsylvan_qasm qsylvan popt)

# allow the parser to be included in other places as a library
# (without compiling eval_expr.cpp again)
add_library(qsylvan_qasm_parser)
target_sources(qsylvan_qasm_parser 
    PRIVATE
        simple_parser.cpp
        parse_math/eval_expr.cpp
    PUBLIC 
        simple_parser.h)

add_executable(run_qasm_on_qmdd run_qasm_on_qmdd.c)
target_link_libraries(run_qasm_on_qmdd qsylvan qsylvan_qasm_parser)

