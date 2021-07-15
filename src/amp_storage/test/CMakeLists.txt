
add_executable(test_amp_storage test_amp_storage.c)
target_link_libraries(test_amp_storage amp_storage)

#add_executable(test_mpfr_map test_mpfr_map.c)
#target_link_libraries(test_mpfr_map amp_storage)

add_test(test_amp_storage test_amp_storage)
#add_test(test_mpfr_map test_mpfr_map)