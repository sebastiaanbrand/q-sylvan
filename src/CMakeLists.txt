add_library(qsylvan STATIC)
add_library(qsylvan::sylvan ALIAS qsylvan)

target_sources(qsylvan
  PRIVATE
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
    sylvan_hash.c
    sylvan_ldd.c
    sylvan_mt.c
    sylvan_mtbdd.c
    sylvan_obj.cpp
    sylvan_refs.c
    sylvan_sl.c
    sylvan_stats.c
    sylvan_table.c
    sylvan_zdd.c
  PUBLIC
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
    sylvan_hash.h
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
    sylvan_zdd.h
    sylvan_zdd_int.h
)

add_subdirectory(edge_weight_storage)
set_target_properties(qsylvan PROPERTIES VERSION ${sylvan_VERSION} SOVERSION ${sylvan_VERSION_MAJOR})
target_compile_features(qsylvan PUBLIC c_std_11 cxx_std_11)
target_compile_options(qsylvan PRIVATE -Wall -Wextra -Werror -fno-strict-aliasing -Wno-deprecated)
target_include_directories(qsylvan PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(qsylvan PUBLIC m pthread lace edge_weight_storage)

find_package(GMP REQUIRED)
target_sources(qsylvan PRIVATE sylvan_gmp.c PUBLIC sylvan_gmp.h)
target_include_directories(qsylvan PRIVATE ${GMP_INCLUDE_DIR})
target_link_libraries(qsylvan PUBLIC ${GMP_LIBRARIES})

find_package(MPC REQUIRED)
target_sources(qsylvan PRIVATE sylvan_mpc.c PUBLIC sylvan_mpc.h)
target_include_directories(qsylvan PRIVATE ${MPC_INCLUDE_DIR})
target_link_libraries(qsylvan PUBLIC ${MPC_LIBRARIES})

# Do we want to collect BDD statistics?
option(SYLVAN_STATS "Let Sylvan collect statistics at runtime" OFF)
if(SYLVAN_STATS)
    set_target_properties(qsylvan PROPERTIES COMPILE_DEFINITIONS "SYLVAN_STATS")
endif()

install(TARGETS qsylvan DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(FILES ${HEADERS} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
