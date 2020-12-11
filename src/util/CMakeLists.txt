

cmake_minimum_required(VERSION 2.7)

include_directories(..)


add_library(util STATIC
        cmap.c cmap.h
        rmap.c rmap.h
        fast_hash.h fast_hash.c
        MurmurHash3.h MurmurHash3.c
        util.h)
