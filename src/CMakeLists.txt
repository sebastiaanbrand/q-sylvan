set(SOURCES
    lace.c
    sha2.c
    qsylvan_gates.c
    qsylvan_simulator.c
    sylvan_aadd.c
    sylvan_bdd.c
    sylvan_cache.c
    sylvan_common.c
    sylvan_edge_weights.c
    sylvan_edge_weights_complex.c
    sylvan_gmp.c
    sylvan_ldd.c
    sylvan_mt.c
    sylvan_mtbdd.c
    sylvan_obj.cpp
    sylvan_refs.c
    sylvan_sl.c
    sylvan_stats.c
    sylvan_table.c
)

set(HEADERS
    lace.h
    sylvan.h
    qsylvan.h
    qsylvan_gates.h
    qsylvan_simulator.h
    sylvan_aadd.h
    sylvan_aadd_int.h
    sylvan_bdd.h
    sylvan_cache.h
    sylvan_config.h
    sylvan_common.h
    sylvan_edge_weights.h
    sylvan_edge_weights_complex.h
    sylvan_gmp.h
    sylvan_int.h
    sylvan_ldd.h
    sylvan_ldd_int.h
    sylvan_mt.h
    sylvan_mtbdd.h
    sylvan_mtbdd_int.h
    sylvan_obj.hpp
    sylvan_stats.h
    sylvan_table.h
    sylvan_tls.h
)

add_subdirectory(edge_weight_storage)

option(BUILD_SHARED_LIBS "Enable/disable creation of shared libraries" ON)
option(BUILD_STATIC_LIBS "Enable/disable creation of static libraries" ON)

add_library(qsylvan ${SOURCES})

find_package(GMP REQUIRED)

include_directories(${GMP_INCLUDE_DIR})

target_link_libraries(qsylvan m pthread ${GMP_LIBRARIES} edge_weight_storage)

if(UNIX AND NOT APPLE)
    target_link_libraries(qsylvan rt)
endif()

option(SYLVAN_STATS "Collect statistics" OFF)
if(SYLVAN_STATS)
    set_target_properties(qsylvan PROPERTIES COMPILE_DEFINITIONS "SYLVAN_STATS")
endif()

install(TARGETS qsylvan DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(FILES ${HEADERS} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
