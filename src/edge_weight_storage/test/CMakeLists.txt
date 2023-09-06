
add_executable(test_wgt_storage test_wgt_storage.c)
target_link_libraries(test_wgt_storage edge_weight_storage m) # clang-14 wants m (math lib) to be linked explicitly

#add_executable(test_mpfr_map test_mpfr_map.c)
#target_link_libraries(test_mpfr_map edge_weight_storage)

add_test(test_wgt_storage test_wgt_storage)
#add_test(test_mpfr_map test_mpfr_map)