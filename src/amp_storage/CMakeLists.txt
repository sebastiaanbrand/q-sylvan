

cmake_minimum_required(VERSION 2.7)

include_directories(..)


add_library(amp_storage SHARED
        amp_storage_interface.c amp_storage_interface.h
        tree_map.cpp tree_map.h
        mpfr_tree_map.cpp mpfr_tree_map.h
        mpreal_tree_map.cpp mpreal_tree_map.h
        cmap.c cmap.h
        rmap.c rmap.h
        fast_hash.h fast_hash.c
        flt.h
        MurmurHash3.h MurmurHash3.c
        util.h)

target_link_libraries(amp_storage "-lquadmath" mpfr)

