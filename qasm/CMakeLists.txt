# 1. Compile qasm parser as separate library
#   (to allow the parser to be included in other places as a library, without 
#   compiling eval_expr.cpp again)
add_library(qsylvan_qasm_parser)
target_sources(qsylvan_qasm_parser 
    PRIVATE
        qsylvan_qasm_parser.cpp
        parse_math/eval_expr.cpp
    PUBLIC 
        qsylvan_qasm_parser.h)

# 2. QMDD QASM simulator
add_executable(run_qasm_on_qmdd run_qasm_on_qmdd.c)
target_link_libraries(run_qasm_on_qmdd qsylvan qsylvan_qasm_parser)

# 3. MTBDD QASM simulator
add_executable(run_qasm_on_mtbdd run_qasm_on_mtbdd.c)
target_link_libraries(run_qasm_on_mtbdd qsylvan qsylvan_qasm_parser)
