add_executable(test_basic test_basic.c)
target_link_libraries(test_basic sylvan)

add_executable(test_cxx test_cxx.cpp)
target_link_libraries(test_cxx sylvan stdc++)

add_executable(qdd_unit_tests qdd_unit_tests.c)
target_link_libraries(qdd_unit_tests sylvan)

add_executable(matrix_qdd_unit_tests matrix_qdd_unit_tests.c)
target_link_libraries(matrix_qdd_unit_tests sylvan)

add_executable(bench_qdd bench_qdd.c)
target_link_libraries(bench_qdd sylvan)

include(CheckIncludeFiles)
check_include_files("gperftools/profiler.h" HAVE_PROFILER)

if(HAVE_PROFILER)
    set_target_properties(bench_qdd PROPERTIES COMPILE_DEFINITIONS "HAVE_PROFILER")
    target_link_libraries(bench_qdd profiler)
endif()

add_test(test_cxx test_cxx)
add_test(test_basic test_basic)
