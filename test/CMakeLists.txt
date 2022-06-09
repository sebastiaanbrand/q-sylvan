# sylvan basic test
add_executable(test_basic test_basic.c)
target_link_libraries(test_basic qsylvan)
add_test(test_basic test_basic)

# sylvan c++ interface test
add_executable(test_cxx test_cxx.cpp)
target_link_libraries(test_cxx qsylvan stdc++)
add_test(test_cxx test_cxx)

# test_qdd_normalization
add_executable(test_qdd_normalization test_qdd_normalization.c)
target_link_libraries(test_qdd_normalization qsylvan)
add_test(test_qdd_normalization test_qdd_normalization)

# test_qdd_basics
add_executable(test_qdd_basics test_qdd_basics.c)
target_link_libraries(test_qdd_basics qsylvan)
add_test(test_qdd_basics test_qdd_basics)

# test_qdd_gates
add_executable(test_qdd_gates test_qdd_gates.c)
target_link_libraries(test_qdd_gates qsylvan)
add_test(test_qdd_gates test_qdd_gates)

# test_qdd_circuits
add_executable(test_qdd_circuits test_qdd_circuits.c)
target_link_libraries(test_qdd_circuits qsylvan)
add_test(test_qdd_circuits test_qdd_circuits)

# test_qdd_matrix
add_executable(test_qdd_matrix test_qdd_matrix.c)
target_link_libraries(test_qdd_matrix qsylvan)
add_test(test_qdd_matrix test_qdd_matrix)

