cmake_minimum_required(VERSION 3.1)

include_directories(..)

add_library(edge_weight_storage SHARED
        wgt_storage_interface.c wgt_storage_interface.h
        cmap.c cmap.h
        fast_hash.h fast_hash.c
        flt.h
        MurmurHash3.h MurmurHash3.c
        util.h)


option(SYLVAN_QUADMATH "128 bit arithmetic" OFF)
if(SYLVAN_QUADMATH)
    target_link_libraries(edge_weight_storage "-lquadmath")
endif()


add_subdirectory(test)
