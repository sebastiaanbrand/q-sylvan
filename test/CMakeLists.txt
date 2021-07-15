# sylvan basic test
add_executable(test_basic test_basic.c)
target_link_libraries(test_basic sylvan)
add_test(test_basic test_basic)

# sylvan c++ interface test
add_executable(test_cxx test_cxx.cpp)
target_link_libraries(test_cxx sylvan stdc++)
add_test(test_cxx test_cxx)

# test_qdd_basics
add_executable(test_qdd_basics test_qdd_basics.c)
target_link_libraries(test_qdd_basics sylvan)
add_test(test_qdd_basics test_qdd_basics)

# test_qdd_gates
add_executable(test_qdd_gates test_qdd_gates.c)
target_link_libraries(test_qdd_gates sylvan)
add_test(test_qdd_gates test_qdd_gates)

# test_qdd_circuits
add_executable(test_qdd_circuits test_qdd_circuits.c)
target_link_libraries(test_qdd_circuits sylvan)
add_test(test_qdd_circuits test_qdd_circuits)

# test_qdd_matrix
add_executable(test_qdd_matrix test_qdd_matrix.c)
target_link_libraries(test_qdd_matrix sylvan)
add_test(test_qdd_matrix test_qdd_matrix)

