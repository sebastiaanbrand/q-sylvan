

cmake_minimum_required(VERSION 2.7)

include_directories(..)


add_library(util #STATIC
        tree_map.cpp tree_map.h
        cmap.c cmap.h
        rmap.c rmap.h
        fast_hash.h fast_hash.c
        flt.h
        MurmurHash3.h MurmurHash3.c
        util.h)

target_link_libraries(util "-lquadmath")
