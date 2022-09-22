# sylvan basic test
add_executable(test_basic test_basic.c)
target_link_libraries(test_basic qsylvan)
add_test(test_basic test_basic)

# sylvan c++ interface test
add_executable(test_cxx test_cxx.cpp)
target_link_libraries(test_cxx qsylvan stdc++)
add_test(test_cxx test_cxx)

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

