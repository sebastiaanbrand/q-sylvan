cmake_minimum_required(VERSION 3.1)

include_directories(..)

add_library(edge_weight_storage SHARED
        wgt_storage_interface.c wgt_storage_interface.h
        tree_map.cpp tree_map.h
        cmap.c cmap.h
        rmap.c rmap.h
        fast_hash.h fast_hash.c
        flt.h
        MurmurHash3.h MurmurHash3.c
        util.h)

target_link_libraries(edge_weight_storage "-lquadmath")

add_subdirectory(test)
