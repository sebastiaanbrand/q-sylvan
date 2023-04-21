include(CTest)

# Build test_basic, test_cxx and test_zdd executables
add_executable(test_basic)
target_sources(test_basic PRIVATE test_basic.c)
target_link_libraries(test_basic PRIVATE qsylvan::sylvan)
target_compile_features(test_basic PRIVATE c_std_11)
target_compile_options(test_basic PRIVATE -Wall -Wextra -Werror -Wno-deprecated)

add_executable(test_mtbdd)
target_sources(test_mtbdd PRIVATE test_mtbdd.c)
target_link_libraries(test_mtbdd PRIVATE qsylvan::sylvan)
target_compile_features(test_mtbdd PRIVATE c_std_11)
target_compile_options(test_mtbdd PRIVATE -Wall -Wextra -Werror -Wno-deprecated)

add_executable(test_cxx)
target_sources(test_cxx PRIVATE test_cxx.cpp)
target_link_libraries(test_cxx PRIVATE qsylvan::sylvan)
target_compile_features(test_cxx PRIVATE cxx_std_11)
target_compile_options(test_cxx PRIVATE -Wall -Wextra -Werror -Wno-deprecated)

add_executable(test_zdd)
target_sources(test_zdd PRIVATE test_zdd.c)
target_link_libraries(test_zdd PRIVATE qsylvan::sylvan)
target_compile_features(test_zdd PRIVATE c_std_11)
target_compile_options(test_zdd PRIVATE -Wall -Wextra -Werror -Wno-deprecated)

# tests for mtb, bdd and ldd 
add_test(test_basic test_basic)
add_test(test_mtbdd test_mtbdd)
add_test(test_cxx test_cxx)
add_test(test_zdd test_zdd)


# test_qmdd_normalization
add_executable(test_qmdd_normalization test_qmdd_normalization.c)
target_link_libraries(test_qmdd_normalization qsylvan)
add_test(test_qmdd_normalization test_qmdd_normalization)

# test_qmdd_basics
add_executable(test_qmdd_basics test_qmdd_basics.c)
target_link_libraries(test_qmdd_basics qsylvan)
add_test(test_qmdd_basics test_qmdd_basics)

# test_qmdd_gates
add_executable(test_qmdd_gates test_qmdd_gates.c)
target_link_libraries(test_qmdd_gates qsylvan)
add_test(test_qmdd_gates test_qmdd_gates)

# test_qmdd_circuits
add_executable(test_qmdd_circuits test_qmdd_circuits.c)
target_link_libraries(test_qmdd_circuits qsylvan)
add_test(test_qmdd_circuits test_qmdd_circuits)

# test_qmdd_matrix
add_executable(test_qmdd_matrix test_qmdd_matrix.c)
target_link_libraries(test_qmdd_matrix qsylvan)
add_test(test_qmdd_matrix test_qmdd_matrix)

# test_qmdd_gc
add_executable(test_qmdd_gc test_qmdd_gc.c ../examples/grover.c)
target_link_libraries(test_qmdd_gc qsylvan)
add_test(test_qmdd_gc test_qmdd_gc)

