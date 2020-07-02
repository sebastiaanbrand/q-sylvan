add_executable(test_basic test_basic.c)
target_link_libraries(test_basic sylvan)

add_executable(test_cxx test_cxx.cpp)
target_link_libraries(test_cxx sylvan stdc++)

add_executable(qdd_unit_tests qdd_unit_tests.c)
target_link_libraries(qdd_unit_tests sylvan)

add_executable(bench_qdd bench_qdd.c)
target_link_libraries(bench_qdd sylvan)

add_test(test_cxx test_cxx)
add_test(test_basic test_basic)
