add_executable(test_basic test_basic.c)
target_link_libraries(test_basic sylvan)

add_executable(test_cxx test_cxx.cpp)
target_link_libraries(test_cxx sylvan stdc++)

add_executable(qdd_unit_tests qdd_unit_tests.c)
target_link_libraries(qdd_unit_tests sylvan)

# TODO: move tests for specific tables to the util folder
add_executable(mpreal_unit_tests mpreal_unit_tests.cpp)
target_link_libraries(mpreal_unit_tests util)

add_executable(matrix_qdd_unit_tests matrix_qdd_unit_tests.c)
target_link_libraries(matrix_qdd_unit_tests sylvan)

add_test(test_cxx test_cxx)
add_test(test_basic test_basic)
