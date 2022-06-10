cmake_minimum_required(VERSION 3.1)

project(qsylvan VERSION 1.0.0) # Based on Sylvan 1.4.0

set(PROJECT_DESCRIPTION "Q-Sylvan, QMDD extension of Sylvan")
set(PROJECT_URL "https://github.com/sebastiaanbrand/q-sylvan")

message(STATUS "CMake build configuration for Q-Sylvan ${PROJECT_VERSION}")
enable_language(C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Some info
message(STATUS "OS: ${CMAKE_SYSTEM_NAME}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

set(CMAKE_C_FLAGS "-O2 -Wextra -Wall -fno-strict-aliasing ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-O2 -Wextra -Wall -fno-strict-aliasing -Wno-deprecated ${CMAKE_CXX_FLAGS}")


option(WITH_COVERAGE "Add generation of test coverage" ON)
if(WITH_COVERAGE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -coverage")

    find_program(GCOV_PATH gcov)
    find_program(LCOV_PATH lcov)
    find_program(GENHTML_PATH genhtml)

    add_custom_target(coverage
        # Cleanup lcov
        ${LCOV_PATH} --directory . --zerocounters
        # Run tests
        COMMAND make test
        # Capture counters
        COMMAND ${LCOV_PATH} --gcov-tool ${GCOV_PATH} --directory . --capture --output-file coverage.info
        COMMAND ${LCOV_PATH} --remove coverage.info 'test/*' '/usr/*' 'examples/*' 'src/lace*' 'src/sha2.c' --output-file coverage.info.cleaned
        COMMAND ${GENHTML_PATH} -o coverage coverage.info.cleaned
        COMMAND ${CMAKE_COMMAND} -E remove coverage.info coverage.info.cleaned
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endif()

include_directories(${CMAKE_SOURCE_DIR}/src)
include(GNUInstallDirs)
include(CTest)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(qasm)
add_subdirectory(examples)

set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "QMDD extension of Sylvan")
set(CPACK_PACKAGE_VENDOR "Sebastiaan Brand")
set(CPACK_PACKAGE_CONTACT "Sebastiaan Brand <s.o.brand@liacs.leidenuniv.nl>")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_PATCH_VERSION})
include(CPack)

# pkg-config
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/qsylvan.pc.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/qsylvan.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qsylvan.pc" DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
