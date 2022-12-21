set(SOURCES
    grover_cnf.c
    grover.c
    random_circuit.c
    shor.c
    supremacy.c
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # add argp library for OSX
    find_package(Argp REQUIRED)
    include_directories(${ARGP_INCLUDE_PATH})
endif()

macro(set_compilation_settings NAME)
    set_target_properties(${NAME} PROPERTIES COMPILE_OPTIONS "-Wno-unused")
    target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR})
    target_link_libraries(${NAME} PRIVATE qsylvan::sylvan ${ARGP_LIBRARIES})
    target_compile_features(${NAME} PRIVATE c_std_11 cxx_std_11)
    target_compile_options(${NAME} PRIVATE -Wall -Wextra -Werror -Wno-deprecated)
endmacro(set_compilation_settings)

macro(add_example NAME SOURCE)
    add_executable(${NAME} ${SOURCE})
    set_compilation_settings(${NAME})
endmacro(add_example)

add_example(bell_state bell_state.c)
# target_sources(bell_state PRIVATE getrss.c getrss.h)

add_example(alg_run alg_run.c)
target_sources(alg_run PRIVATE ${SOURCES})

add_example(test_algs test_algs.c)
target_sources(test_algs PRIVATE ${SOURCES})
